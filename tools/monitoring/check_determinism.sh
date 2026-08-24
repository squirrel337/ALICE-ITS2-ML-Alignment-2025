#!/bin/bash
# Run the same job twice from two pristine copies and compare.
#
# The module is not bit-reproducible run to run: two byte-identical job trees,
# same machine, same binary, produce different numbers from the very first
# computed quantity onward -- the epoch -1 cost, which is a pure forward
# evaluation with every sensor weight at zero. This script measures the size of
# that spread so it can be used as the acceptance band when comparing two
# versions of the module, since bit-identity is not achievable.
#
#   check_determinism.sh <composed-job-dir> [n]
#
# Copies the job n times (default 2), runs each once in its own directory, and
# reports the first epoch -1 cost from each plus where the logs first diverge.
# Use a small nDATA -- the spread is visible immediately and a full job is hours.
set -u
SRC=${1:?usage: check_determinism.sh <job-dir> [n]}
N=${2:-2}
BASE=$(cd "$(dirname "$SRC")" && pwd)/$(basename "$SRC")
WORK=${BASE}.determinism
rm -rf "$WORK"; mkdir -p "$WORK"

echo "copying $N pristine job trees ..."
for i in $(seq 1 "$N"); do
  cp -a "$BASE" "$WORK/r$i"
  # a run consumes its own directory (it mv's UpdateSensorsList.txt and
  # TrendingNetwork/ into MLPTrain/), so every repetition needs its own copy --
  # re-running in place is not the same experiment.
  rm -rf "$WORK/r$i/MLPTrain" "$WORK/r$i"/MLPTrain_Step[0-9]*/ 2>/dev/null
done

for i in $(seq 1 "$N"); do
  ( cd "$WORK/r$i" && root -l -b < batch_train > train.log 2>&1 )
  printf "r%-3s %s\n" "$i" "$(grep -m1 -oE 'Fit \+ CHSYM = [0-9.]+ \+ [0-9.]+' "$WORK/r$i/train.log")"
done

echo
echo "first divergence between r1 and r2:"
sed "s|/r1/|/JOB/|g" "$WORK/r1/train.log" > "$WORK/a.log"
sed "s|/r2/|/JOB/|g" "$WORK/r2/train.log" > "$WORK/b.log"
if diff -q "$WORK/a.log" "$WORK/b.log" >/dev/null; then
  echo "  none -- the two runs are identical"
else
  diff "$WORK/a.log" "$WORK/b.log" | head -6
fi
echo
echo "weights, if the runs got that far:"
for i in 1 2; do
  w=$(ls "$WORK/r$i"/MLPTrain_Step*/weights/weights_Epoch_At_0.txt 2>/dev/null | head -1)
  [ -n "$w" ] && echo "  r$i $w"
done
echo "compare them with tools/monitoring/compare_weight_files.py"
