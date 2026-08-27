#!/bin/bash
# node 252: REDO magnetic L=12, U=2, delta 0.5/0.6/0.7 (lost to the maintenance reboot).
# Driver now writes after EACH delta, so a reboot costs at most one delta rather than the block.
S=~/Desktop/checkerboard-altermagnet/code
ssh 252 'mkdir -p ~/mag_hi'
scp "$S/cpqmc.py" "$S/checkerboard.py" "$S/unified_scan.py" "$S/checkerboard_mag_scan.py" 252:~/mag_hi/
ssh -f 252 'cd ~/mag_hi && SIZES=12 US=2 DELTAS=0.5,0.6,0.7 NPROC=32 NSEED=6 nohup /opt/anaconda3/bin/python checkerboard_mag_scan.py </dev/null > mag_hi.log 2>&1 &'
sleep 1; echo "== 252: magnetic L=12 U=2, delta 0.5/0.6/0.7 (incremental writes) =="
