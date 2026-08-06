#!/bin/bash
# node 251: PAIRING cube, L=12, U=8, HIGH delta = 0.5, 0.6, 0.7
# Separate dir ~/pair_hi/ so it cannot overwrite the existing delta<=0.4 pairing files.
S=~/Desktop/checkerboard-altermagnet/code
ssh 251 'mkdir -p ~/pair_hi'
scp "$S/cpqmc.py" "$S/checkerboard.py" "$S/checkerboard_fss_full.py" 251:~/pair_hi/
ssh -f 251 'cd ~/pair_hi && SIZES=12 US=8 DELTAS=0.5,0.6,0.7 NPROC=32 NSEED=6 nohup ~/anaconda3/bin/python checkerboard_fss_full.py </dev/null > pair_hi.log 2>&1 &'
sleep 1; echo "== 251: pairing L=12, U=8, delta 0.5/0.6/0.7 =="
