#!/bin/bash
# Repeat one composed job N times from N pristine copies and record the numbers
# that vary between runs.
#
#   run_to_run.sh <composed-job-dir> <label> [n] [parallel]
#
# A run consumes its own directory -- run_train_circle.C moves UpdateSensorsList.txt
# and TrendingNetwork/ into MLPTrain/ -- so every repetition gets its own copy.
# Re-running in place is a different experiment and gives a different answer.
#
# Writes <job-dir>.rtr-<label>/r<i>/train.log and a summary at
# <job-dir>.rtr-<label>/costs.tsv with one line per run:
#
#   <i>  <fit>  <chsym>  <status>
#
# Intended for nEPOCH=0, which evaluates without training: the epoch -1 cost is
# a pure forward pass with every sensor weight at zero, so anything that differs
# between two runs of it is non-determinism and nothing else.
set -u

SRC=${1:?usage: run_to_run.sh <job-dir> <label> [n] [parallel]}
LABEL=${2:?usage: run_to_run.sh <job-dir> <label> [n] [parallel]}
N=${3:-30}
# Parallelism is bounded by memory, not by cores. One run cling-JITs the whole
# module and then holds the event sample: measured 7.8 GB resident at nDATA=4000,
# still climbing through the epoch. Two of those do not fit under a 13 GB cgroup.
#
# Overcommitting does not fail loudly. The run is OOM-killed at its peak, which is
# just after the network is built, and the log stops at exactly the line a run with
# nEPOCH=0 stops at -- no error, no cost, no weights, in about a third of the time.
# Eight repetitions were lost to this before the cgroup limit was checked
# (/sys/fs/cgroup/memory/.../memory.limit_in_bytes, and memory.max_usage_in_bytes
# pinned to it). Check the budget before raising this.
PAR=${4:-1}

BASE=$(cd "$(dirname "$SRC")" && pwd)/$(basename "$SRC")
WORK="${BASE}.rtr-${LABEL}"

# The job reads its starting point from MLPTrain_Step<FIRST_STEP-1>: run_train_circle.C
# calls SetPrevUSL, SetPrevWeight and SetPrevWeightDU on that directory. It must
# survive into every copy. Only the *output* of a previous run may be cleared --
# MLPTrain/ (a crashed leftover) and MLPTrain_Step<n> for n >= FIRST_STEP, which
# line 74 produces by renaming MLPTrain at the end.
#
# Deleting the seed does not fail loudly: the job builds the whole network, then
# stops silently right after "[ITS OB2]" with a 25714-line log, no error, no cost
# and no weights, in about half the time a real run takes.
FIRST=$(awk -F= '/^first=/{print $2; exit}' "$BASE/run_steps.sh" 2>/dev/null)
[ -n "${FIRST:-}" ] || { echo "cannot read first= from $BASE/run_steps.sh" >&2; exit 1; }
SEED="MLPTrain_Step$((FIRST - 1))"
[ -d "$BASE/$SEED" ] || {
  echo "$BASE has no $SEED -- the job has no starting point; recompose with PARAMS_ARCHIVE set" >&2
  exit 1; }

rm -rf "$WORK"; mkdir -p "$WORK"
echo "[$LABEL] copying $N pristine job trees into $WORK  (seed $SEED kept)"
for i in $(seq 1 "$N"); do
  cp -a "$BASE" "$WORK/r$i"
  rm -rf "$WORK/r$i/MLPTrain"
  for d in "$WORK/r$i"/MLPTrain_Step[0-9]*; do
    [ -d "$d" ] || continue
    n=${d##*/MLPTrain_Step}
    [ "$n" -ge "$FIRST" ] 2>/dev/null && rm -rf "$d"
  done
done

one() {   # index
  local i=$1 d="$WORK/r$1" s e
  s=$(date +%s)
  ( cd "$d" && root -l -b < batch_train > train.log 2>&1 )
  e=$(date +%s)
  echo "$((e - s))" > "$d/elapsed"
  echo "[$LABEL] r$i done in $((e - s))s"
}

echo "[$LABEL] running $N repetitions, $PAR at a time"
running=0
for i in $(seq 1 "$N"); do
  one "$i" &
  running=$((running + 1))
  if [ "$running" -ge "$PAR" ]; then wait -n 2>/dev/null || wait; running=$((running - 1)); fi
done
wait

# What a repetition is worth is its epoch-0 weight file, not the cost line.
#
# nEPOCH must be at least 1. The epoch loop at YMultiLayerPerceptron.cxx:1748 is
# "for (iepoch = 0; iepoch < nEpoch; ...)", so nEPOCH=0 never enters it: the job
# builds the network, prints nothing, and exits. nEPOCH=1 runs epoch 0 and writes
# weights_Epoch_At_0.txt, which is the parameter set to compare.
#
# Both the cost and the weights are usable, and both depend on the seed being
# there. A run whose MLPTrain_Step<FIRST_STEP-1> is missing prints EPOCH-1 as -nan
# and EPOCH0 as 0 and still writes a weight file, so keying success on the cost
# alone accepts a corrupted run as a measurement. The weight file is the criterion
# because it is what the job exists to produce; the cost is recorded alongside it.
: > "$WORK/costs.tsv"
for i in $(seq 1 "$N"); do
  line=$(grep -m1 -oE 'COSTMONITOR\[TRAINING\] EPOCH0 Fit \+ CHSYM = [0-9.eE+-]+ \+ [0-9.eE+-]+' \
         "$WORK/r$i/train.log" 2>/dev/null)
  fit=$(echo "$line" | awk '{print $7}')
  chs=$(echo "$line" | awk '{print $9}')
  w=$(ls "$WORK/r$i"/MLPTrain_Step*/weights/weights_Epoch_At_0.txt 2>/dev/null | head -1)
  st=ok; [ -z "$w" ] && st=noweights
  printf '%s\t%s\t%s\t%s\t%s\n' "$i" "${fit:-NA}" "${chs:-NA}" "$st" "${w:-NA}" >> "$WORK/costs.tsv"
done

ok=$(awk -F'\t' '$4=="ok"' "$WORK/costs.tsv" | wc -l)
echo "[$LABEL] $ok/$N repetitions produced epoch-0 weights -> $WORK/costs.tsv"
[ "$ok" -lt "$N" ] && awk -F'\t' '$4!="ok"{print "  r"$1" produced no weights"}' "$WORK/costs.tsv"
echo "[$LABEL] analyse with: tools/monitoring/weight_spread.py out.png $LABEL=$WORK/costs.tsv"
