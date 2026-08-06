#!/bin/bash
# node 250 (32 physical cores) -> all sizes L=8,10,12 at U=0
S=~/Desktop/checkerboard-altermagnet/code
ssh 250 'mkdir -p ~/cb_run'
scp "$S/cpqmc.py" "$S/checkerboard.py" "$S/checkerboard_fss_full.py" 250:~/cb_run/
ssh -f 250 'cd ~/cb_run && SIZES=8,10,12 US=0 NPROC=32 NSEED=6 nohup ~/anaconda3/bin/python checkerboard_fss_full.py </dev/null > fss_full.log 2>&1 &'
sleep 1; echo "== 250 launched: L=8,10,12 at U=0 =="
