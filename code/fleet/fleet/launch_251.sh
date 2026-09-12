#!/bin/bash
# node 251 (32 physical cores) -> all sizes L=8,10,12 at U=8
S=~/Desktop/checkerboard-altermagnet/code
ssh 251 'mkdir -p ~/cb_run'
scp "$S/cpqmc.py" "$S/checkerboard.py" "$S/checkerboard_fss_full.py" 251:~/cb_run/
ssh -f 251 'cd ~/cb_run && SIZES=8,10,12 US=8 NPROC=32 NSEED=6 nohup ~/anaconda3/bin/python checkerboard_fss_full.py </dev/null > fss_full.log 2>&1 &'
sleep 1; echo "== 251 launched: L=8,10,12 at U=8 =="
