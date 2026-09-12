#!/bin/bash
# node 255: magnetic scan, L=12 ONLY (never-trust-L6 rule), U=0,2,4
S=~/Desktop/checkerboard-altermagnet/code
P=~/Desktop/Susceptibility/qmc-platform-master/pyqmc
ssh 255 'mkdir -p ~/mag_run'
scp "$S/checkerboard.py" "$S/checkerboard_mag_scan.py" 255:~/mag_run/
scp "$P/cpqmc.py" "$P/unified_scan.py" 255:~/mag_run/
ssh -f 255 'cd ~/mag_run && SIZES=12 US=0,2,4 NPROC=32 NSEED=6 nohup /opt/anaconda3/bin/python checkerboard_mag_scan.py </dev/null > mag.log 2>&1 &'
sleep 1; echo "== 255 launched: magnetic L=12, U=0,2,4 =="
