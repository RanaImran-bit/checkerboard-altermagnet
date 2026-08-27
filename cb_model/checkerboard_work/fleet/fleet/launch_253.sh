#!/bin/bash
# node 253 (32 physical cores) -> all sizes L=8,10,12 at U=2
S=~/Desktop/checkerboard-altermagnet/code
ssh 253 'mkdir -p ~/cb_run'
scp "$S/cpqmc.py" "$S/checkerboard.py" "$S/checkerboard_fss_full.py" 253:~/cb_run/
ssh -f 253 'cd ~/cb_run && SIZES=8,10,12 US=2 NPROC=32 NSEED=6 nohup /opt/anaconda3/bin/python checkerboard_fss_full.py </dev/null > fss_full.log 2>&1 &'
sleep 1; echo "== 253 launched: L=8,10,12 at U=2 =="
