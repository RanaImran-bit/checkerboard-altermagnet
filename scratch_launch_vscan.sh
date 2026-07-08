#!/bin/bash
# V-SCAN: does the d-wave vertex rise/kink/fall when S_AM switches on?
# Fixed manuscript bands (t1=1,t2=1.75,t3=0.85,t4=0.65) + uxy=1; scan neighbour v=0..1.
# Sharded round-robin. Usage: launch_vscan.sh SHARD NSHARD
set -u
SHARD=${1:-0}; NSHARD=${2:-1}
cd ~/qmc_twoorb
PY=""; for c in ~/miniconda3/bin/python3 /usr/bin/python3 python3; do "$c" -c "import numpy" 2>/dev/null && { PY="$c"; break; }; done
[ -z "$PY" ] && { echo NO_NUMPY; exit 2; }
OUT=~/qmc_twoorb/out_v; mkdir -p $OUT
export OMP_NUM_THREADS=1 KMP_DUPLICATE_LIB_OK=TRUE
T1=1; T2=1.75; T3=0.85; T4=0.65; UXY=1
VLIST=(0.0 0.2 0.4 0.6 0.8 1.0)
idx=0; launched=0
launch(){ local tag=$1; shift; setsid nohup "$@" >"$OUT/$tag.log" 2>&1 </dev/null & }

for vi in 0 1 2 3 4 5; do
  V=${VLIST[$vi]}
  # (1) DQMC 6x6, U=4, beta=4 (exact; sign marginal ~0.5, tracked)
  if [ $((idx%NSHARD)) -eq $SHARD ]; then
    launch "v_dqmc_L6U4_v${vi}" $PY code/dqmc_py/two_orb_ft_scan.py --lx 6 --ly 6 --U 4 --uxy $UXY --v $V \
      --beta 4 --dt 0.1 --t1 $T1 --t3 $T3 --t4 $T4 --method dqmc --nwarm 50 --nmeas 260 --seed 1 --t2list $T2 \
      -o "$OUT/v_dqmc_L6U4_v${vi}.csv"; launched=$((launched+1)); fi
  idx=$((idx+1))
  # (2) DQMC 6x6, U=2, beta=5 (exact, sign~1 clean reference)
  if [ $((idx%NSHARD)) -eq $SHARD ]; then
    launch "v_dqmc_L6U2_v${vi}" $PY code/dqmc_py/two_orb_ft_scan.py --lx 6 --ly 6 --U 2 --uxy $UXY --v $V \
      --beta 5 --dt 0.125 --t1 $T1 --t3 $T3 --t4 $T4 --method dqmc --nwarm 50 --nmeas 260 --seed 1 --t2list $T2 \
      -o "$OUT/v_dqmc_L6U2_v${vi}.csv"; launched=$((launched+1)); fi
  idx=$((idx+1))
  # (3) DQMC 8x8, U=2, beta=4 (size check at sign-safe U)
  if [ $((idx%NSHARD)) -eq $SHARD ]; then
    launch "v_dqmc_L8U2_v${vi}" $PY code/dqmc_py/two_orb_ft_scan.py --lx 8 --ly 8 --U 2 --uxy $UXY --v $V \
      --beta 4 --dt 0.125 --t1 $T1 --t3 $T3 --t4 $T4 --method dqmc --nwarm 40 --nmeas 170 --seed 1 --t2list $T2 \
      -o "$OUT/v_dqmc_L8U2_v${vi}.csv"; launched=$((launched+1)); fi
  idx=$((idx+1))
  # (4) CPQMC 6x6, U=4, T=0 (reaches where the sign dies; cheap settings, full v cost)
  if [ $((idx%NSHARD)) -eq $SHARD ]; then
    launch "v_cpqmc_L6U4_v${vi}" $PY pyqmc/two_orb_uscan.py --lx 6 --ly 6 --nup 36 --ndn 36 \
      --t1 $T1 --t2 $T2 --t3 $T3 --t4 $T4 --uxy $UXY --v $V --dt 0.05 --bp 10 --nw 48 --nequil 50 --nblocks 10 \
      --seed 1 --Ulist 4 -o "$OUT/v_cpqmc_L6U4_v${vi}.csv"; launched=$((launched+1)); fi
  idx=$((idx+1))
  # (5) CP-DQMC 6x6, U=4, beta=4 (constrained sign=1 cross-check; very cheap; even vi only)
  if [ $((vi%2)) -eq 0 ] && [ $((idx%NSHARD)) -eq $SHARD ]; then
    launch "v_cpdqmc_L6U4_v${vi}" $PY code/dqmc_py/two_orb_ft_scan.py --lx 6 --ly 6 --U 4 --uxy $UXY --v $V \
      --beta 4 --dt 0.1 --t1 $T1 --t3 $T3 --t4 $T4 --method cpdqmc --nwarm 15 --nmeas 30 --nw 8 --seed 1 --t2list $T2 \
      -o "$OUT/v_cpdqmc_L6U4_v${vi}.csv"; launched=$((launched+1)); fi
  idx=$((idx+1))
done
sleep 1; echo "VSCAN SHARD $SHARD/$NSHARD PY=$PY launched=$launched procs=$(pgrep -f 'two_orb_(uscan|ft_scan)'|wc -l)"
