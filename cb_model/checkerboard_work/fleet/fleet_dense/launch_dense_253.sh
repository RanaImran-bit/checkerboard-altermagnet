#!/bin/bash
# node 253: dense-filling pairing scan, L=12, U=4, nup=49..60
S=~/Desktop/checkerboard-altermagnet/code
ssh 253 'mkdir -p ~/dense_run'
scp "$S/cpqmc.py" "$S/checkerboard.py" "$S/checkerboard_dense_n.py" 253:~/dense_run/ 2>/dev/null ||   scp ~/Desktop/Susceptibility/qmc-platform-master/pyqmc/cpqmc.py "$S/checkerboard.py" "$S/checkerboard_dense_n.py" 253:~/dense_run/
ssh -f 253 'cd ~/dense_run && L=12 U=4 NUP_MIN=49 NUP_MAX=60 STEP=1 NPROC=32 NSEED=6 nohup /opt/anaconda3/bin/python checkerboard_dense_n.py </dev/null > dense.log 2>&1 &'
sleep 1; echo "== 253 launched: dense n, nup=49..60 (n=$(python3 -c "print(f'{2*49/144:.3f}')")..$(python3 -c "print(f'{2*60/144:.3f}')")) =="
