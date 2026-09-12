#!/bin/bash
# node 252: PAIRING cube, L=12, U=2, HIGH delta = 0.5, 0.6, 0.7
# Separate dir ~/pair_hi/ so it cannot overwrite the existing delta<=0.4 pairing files.
S=~/Desktop/checkerboard-altermagnet/code
ssh 252 'mkdir -p ~/pair_hi'
scp "$S/cpqmc.py" "$S/checkerboard.py" "$S/checkerboard_fss_full.py" 252:~/pair_hi/
ssh -f 252 'cd ~/pair_hi && SIZES=12 US=2 DELTAS=0.5,0.6,0.7 NPROC=32 NSEED=6 nohup /opt/anaconda3/bin/python checkerboard_fss_full.py </dev/null > pair_hi.log 2>&1 &'
sleep 1; echo "== 252: pairing L=12, U=2, delta 0.5/0.6/0.7 =="
