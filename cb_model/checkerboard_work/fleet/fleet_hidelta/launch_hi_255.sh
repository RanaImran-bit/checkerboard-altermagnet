#!/bin/bash
# node 255: magnetic scan, L=12, U=0, HIGH delta = 0.5, 0.6, 0.7
# Separate dir ~/mag_hi/ so it cannot overwrite the existing delta<=0.4 files.
S=~/Desktop/checkerboard-altermagnet/code
ssh 255 'mkdir -p ~/mag_hi'
scp "$S/cpqmc.py" "$S/checkerboard.py" "$S/unified_scan.py" "$S/checkerboard_mag_scan.py" 255:~/mag_hi/
ssh -f 255 'cd ~/mag_hi && SIZES=12 US=0 DELTAS=0.5,0.6,0.7 NPROC=32 NSEED=6 nohup /opt/anaconda3/bin/python checkerboard_mag_scan.py </dev/null > mag_hi.log 2>&1 &'
sleep 1; echo "== 255: magnetic L=12, U=0, delta 0.5/0.6/0.7 =="
