#!/bin/bash
# node 255 (32 physical cores) -> all sizes L=8,10,12 at U=6
S=~/Desktop/checkerboard-altermagnet/code
ssh 255 'mkdir -p ~/cb_run'
scp "$S/cpqmc.py" "$S/checkerboard.py" "$S/checkerboard_fss_full.py" 255:~/cb_run/
ssh -f 255 'cd ~/cb_run && SIZES=8,10,12 US=6 NPROC=32 NSEED=6 nohup /opt/anaconda3/bin/python checkerboard_fss_full.py </dev/null > fss_full.log 2>&1 &'
sleep 1; echo "== 255 launched: L=8,10,12 at U=6 =="
