#!/bin/bash
# node 250: dense-filling pairing scan, L=12, U=4, nup=36..48
S=~/Desktop/checkerboard-altermagnet/code
ssh 250 'mkdir -p ~/dense_run'
scp "$S/cpqmc.py" "$S/checkerboard.py" "$S/checkerboard_dense_n.py" 250:~/dense_run/ 2>/dev/null ||   scp ~/Desktop/Susceptibility/qmc-platform-master/pyqmc/cpqmc.py "$S/checkerboard.py" "$S/checkerboard_dense_n.py" 250:~/dense_run/
ssh -f 250 'cd ~/dense_run && L=12 U=4 NUP_MIN=36 NUP_MAX=48 STEP=1 NPROC=32 NSEED=6 nohup ~/anaconda3/bin/python checkerboard_dense_n.py </dev/null > dense.log 2>&1 &'
sleep 1; echo "== 250 launched: dense n, nup=36..48 (n=$(python3 -c "print(f'{2*36/144:.3f}')")..$(python3 -c "print(f'{2*48/144:.3f}')")) =="
