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
# One run is a cling JIT of the whole module and peaks near 4.6 GB resident, so
# parallelism is bounded by memory, not by cores. Overcommitting gets a run
# OOM-killed silently -- the log simply stops after the network is built and no
# cost line is ever printed, which costs a whole repetition. Budget ~5 GB each.
PAR=${4:-2}

BASE=$(cd "$(dirname "$SRC")" && pwd)/$(basename "$SRC")
WORK="${BASE}.rtr-${LABEL}"

rm -rf "$WORK"; mkdir -p "$WORK"
echo "[$LABEL] copying $N pristine job trees into $WORK"
for i in $(seq 1 "$N"); do
  cp -a "$BASE" "$WORK/r$i"
  rm -rf "$WORK/r$i/MLPTrain" "$WORK/r$i"/MLPTrain_Step[0-9]*/ 2>/dev/null
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

# The first computed quantity of the job: the epoch 0 training cost, printed as
# "COSTMONITOR[TRAINING] EPOCH0 Fit + CHSYM = <fit> + <chsym>". Everything
# downstream inherits its spread.
#
# nEPOCH must be at least 1. The epoch loop at YMultiLayerPerceptron.cxx:1748 is
# "for (iepoch = 0; iepoch < nEpoch; ...)", so nEPOCH=0 never enters it: the job
# builds the network, prints nothing, and exits. nEPOCH=1 runs epoch 0 and writes
# weights_Epoch_At_0.txt, which is the epoch-0 parameter set to compare.
: > "$WORK/costs.tsv"
for i in $(seq 1 "$N"); do
  line=$(grep -m1 -oE 'COSTMONITOR\[TRAINING\] EPOCH0 Fit \+ CHSYM = [0-9.eE+-]+ \+ [0-9.eE+-]+' \
         "$WORK/r$i/train.log" 2>/dev/null)
  fit=$(echo "$line" | awk '{print $7}')
  chs=$(echo "$line" | awk '{print $9}')
  w=$(ls "$WORK/r$i"/MLPTrain_Step*/weights/weights_Epoch_At_0.txt 2>/dev/null | head -1)
  st=ok; [ -z "$fit" ] && st=nocost
  printf '%s\t%s\t%s\t%s\t%s\n' "$i" "${fit:-NA}" "${chs:-NA}" "$st" "${w:-NA}" >> "$WORK/costs.tsv"
done

echo "[$LABEL] summary -> $WORK/costs.tsv"
awk -F'\t' '$4=="ok"{n++; v[n]=$2; if(n==1||$2<lo)lo=$2; if(n==1||$2>hi)hi=$2}
  END{ if(n<2){print "  n="n" (not enough runs)"; exit}
       printf "  n=%d  min=%s  max=%s  range=%.3e  relative=%.3e\n", n, lo, hi, hi-lo, (hi-lo)/lo }' "$WORK/costs.tsv"
