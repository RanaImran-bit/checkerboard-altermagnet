#!/bin/bash
# node 250: magnetic scan, L=12,10,8 (L=6 dropped per never-trust-L6), U=0
# Ships the minimal unified_scan.py SHIM from code/ -- the real pyqmc one pulls in
# agp_dwave and dies with ModuleNotFoundError on the nodes.
S=~/Desktop/checkerboard-altermagnet/code
ssh 250 'mkdir -p ~/mag_run'
scp "$S/cpqmc.py" "$S/checkerboard.py" "$S/unified_scan.py" "$S/checkerboard_mag_scan.py" 250:~/mag_run/
ssh -f 250 'cd ~/mag_run && SIZES=12,10,8 US=0 NPROC=32 NSEED=6 nohup ~/anaconda3/bin/python checkerboard_mag_scan.py </dev/null > mag.log 2>&1 &'
sleep 1; echo "== 250 relaunched: magnetic L=12,10,8 at U=0 =="
