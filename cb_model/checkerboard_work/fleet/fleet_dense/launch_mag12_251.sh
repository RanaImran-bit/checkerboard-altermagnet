#!/bin/bash
# node 251: magnetic scan, L=12 ONLY (never-trust-L6 rule), U=6,8
S=~/Desktop/checkerboard-altermagnet/code
P=~/Desktop/Susceptibility/qmc-platform-master/pyqmc
ssh 251 'mkdir -p ~/mag_run'
scp "$S/checkerboard.py" "$S/checkerboard_mag_scan.py" 251:~/mag_run/
scp "$P/cpqmc.py" "$P/unified_scan.py" 251:~/mag_run/
ssh -f 251 'cd ~/mag_run && SIZES=12 US=6,8 NPROC=32 NSEED=6 nohup ~/anaconda3/bin/python checkerboard_mag_scan.py </dev/null > mag.log 2>&1 &'
sleep 1; echo "== 251 launched: magnetic L=12, U=6,8 =="
