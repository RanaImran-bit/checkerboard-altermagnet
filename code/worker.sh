#!/bin/bash
# worker.sh <joblist> <worker_index> <n_workers> [NP]
# Each worker takes every n_workers-th line, so W workers share one list with no
# lock file and no risk of two workers taking the same job.
set -u
LIST=$1; IDX=$2; NW=$3; NP=${4:-8}
cd "$HOME/Checkerboard_Model" || exit 1
i=0
while read -r L U DEL NF; do
    [ -z "${L:-}" ] && continue
    if [ $((i % NW)) -eq "$IDX" ]; then
        # < /dev/null: without it mpirun/make read the job list from stdin
        bash run_sync.sh "$L" "$U" "$DEL" "$NF" "$NP" 0.3 < /dev/null
    fi
    i=$((i+1))
done < "$LIST"
echo "WORKER $IDX FINISHED"
