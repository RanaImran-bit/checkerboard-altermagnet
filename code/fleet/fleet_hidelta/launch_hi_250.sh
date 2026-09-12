#!/bin/bash
# node 250: magnetic scan, L=12, U=6, HIGH delta = 0.5, 0.6, 0.7
# Separate dir ~/mag_hi/ so it cannot overwrite the existing delta<=0.4 files.
S=~/Desktop/checkerboard-altermagnet/code
ssh 250 'mkdir -p ~/mag_hi'
scp "$S/cpqmc.py" "$S/checkerboard.py" "$S/unified_scan.py" "$S/checkerboard_mag_scan.py" 250:~/mag_hi/
ssh -f 250 'cd ~/mag_hi && SIZES=12 US=6 DELTAS=0.5,0.6,0.7 NPROC=32 NSEED=6 nohup ~/anaconda3/bin/python checkerboard_mag_scan.py </dev/null > mag_hi.log 2>&1 &'
sleep 1; echo "== 250: magnetic L=12, U=6, delta 0.5/0.6/0.7 =="
