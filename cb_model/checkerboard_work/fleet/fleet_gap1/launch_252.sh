#!/bin/bash
# node 252: magnetic L=8, U=6, delta 0-0.4. Fills the finite-size gap so the
# magnetic cube has all five U at L=8 and L=10, matching L=12.
S=~/Desktop/checkerboard-altermagnet/code
ssh 252 'mkdir -p ~/mag_gap1'
scp "$S/cpqmc.py" "$S/checkerboard.py" "$S/unified_scan.py" "$S/checkerboard_mag_scan.py" 252:~/mag_gap1/
ssh -f 252 'cd ~/mag_gap1 && SIZES=8 US=6 DELTAS=0,0.1,0.2,0.3,0.4 NPROC=32 NSEED=6 nohup /opt/anaconda3/bin/python checkerboard_mag_scan.py </dev/null > mag_gap1.log 2>&1 &'
sleep 1; echo "== 252: magnetic L=8 U=6, delta 0-0.4 =="
