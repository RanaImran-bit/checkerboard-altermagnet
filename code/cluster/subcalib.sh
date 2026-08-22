#!/bin/bash
#SBATCH -p all
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=12G
#SBATCH -o logs/%x-%j.out
# Time ONE L=12 chi cell at 1/20 the step count, ON 113new, to measure this
# cluster's per-core speed. The 256 calibration gave 11.3 h for a full cell but
# on an idle machine; the running jobs sit on nodes at CPU_LOAD 119-147 out of 96
# cores. Only a measurement here converts that into an ETA.
#   full:  NEQ=640 + 40*16=640 -> 1280 steps      here: NEQ=32 + 2*16=32 -> 64 (1/20)
B=/opt/data/hpc/ZONES/all/PROJECTS/a1/PRIVATE/jianyu/Imran/checkerboard
mkdir -p "$B/calib" && cd "$B/calib" || exit 1
cp "$B"/pychi/{checkerboard_eqtime.py,checkerboard.py,cpqmc.py} . 2>/dev/null
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export L=12 NUPS=72 US=2 DELTAS=0.3 NSEED=1 NPROC=1
export NW=500 NEQ=32 NBLK=2 BP=16 DT=0.05
s=$(date +%s)
/opt/data/hpc/common/anaconda3/bin/python checkerboard_eqtime.py
e=$(date +%s); d=$((e-s))
echo "ELAPSED_SEC=$d"
echo "SCALED_FULL_HOURS_113NEW=$(awk -v x=$d 'BEGIN{printf "%.2f", x*20/3600}')"
