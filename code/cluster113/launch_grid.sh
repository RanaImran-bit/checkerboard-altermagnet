#!/bin/bash
# launch_grid.sh <L> -- submit the full 7 U x 7 delta half-filling grid.
# Must be run from the checkerboard dir: the SLURM lua filter rejects any
# submission whose real path is not under the PROJECTS zone.
# Partition MUST be `all`; the default `normal` has MaxTime=00:01:00.
L=${1:?usage: launch_grid.sh <L>}
B=/opt/data/hpc/ZONES/all/PROJECTS/a1/PRIVATE/jianyu/Imran/checkerboard
cd "$B" || exit 1
mkdir -p logs
n=0
for U in 0 2 3 3.5 4 4.5 5; do
  for D in 0.1 0.2 0.3 0.4 0.5 0.6 0.7; do
    sbatch -p all -N 1 -n 8 -J cb${L}u${U}d${D} -o logs/%x-%j.out \
           sub113.sh "$L" "$U" "$D" 1.0 8 0.3 >/dev/null 2>&1 && n=$((n+1))
  done
done
echo "submitted $n jobs for L=$L"
