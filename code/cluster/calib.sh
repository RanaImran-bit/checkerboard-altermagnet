#!/bin/bash
# Time ONE L=12 chi cell at the exact 113new settings, with the step count cut 20x,
# then scale. The 113new jobs run 42 tasks on 42 cores, so one task's wall time IS
# the job's wall time -- and there is no other way to bound them, because the
# driver writes only after Pool.map returns and emits no progress whatsoever.
#   full:  NEQ=640 + NBLK*BP = 40*16 = 640  ->  1280 propagation steps
#   here:  NEQ=32  + NBLK*BP =  2*16 =  32  ->    64 propagation steps  (1/20)
cd "$HOME/chi_calib" || exit 1
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export L=12 NUPS=72 US=2 DELTAS=0.3 NSEED=1 NPROC=1
export NW=500 NEQ=32 NBLK=2 BP=16 DT=0.05
s=$(date +%s)
/opt/anaconda3/bin/python checkerboard_eqtime.py > calib.out 2>&1
e=$(date +%s)
d=$((e-s))
echo "ELAPSED_SEC=$d" >> calib.out
echo "SCALED_FULL_HOURS=$(awk -v x=$d 'BEGIN{printf "%.2f", x*20/3600}')" >> calib.out
