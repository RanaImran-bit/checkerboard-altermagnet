#!/bin/bash
# node 254: dense-filling pairing scan, L=12, U=4, nup=61..72
S=~/Desktop/checkerboard-altermagnet/code
ssh 254 'mkdir -p ~/dense_run'
scp "$S/cpqmc.py" "$S/checkerboard.py" "$S/checkerboard_dense_n.py" 254:~/dense_run/ 2>/dev/null ||   scp ~/Desktop/Susceptibility/qmc-platform-master/pyqmc/cpqmc.py "$S/checkerboard.py" "$S/checkerboard_dense_n.py" 254:~/dense_run/
ssh -f 254 'cd ~/dense_run && L=12 U=4 NUP_MIN=61 NUP_MAX=72 STEP=1 NPROC=32 NSEED=6 nohup /opt/anaconda3/bin/python checkerboard_dense_n.py </dev/null > dense.log 2>&1 &'
sleep 1; echo "== 254 launched: dense n, nup=61..72 (n=$(python3 -c "print(f'{2*61/144:.3f}')")..$(python3 -c "print(f'{2*72/144:.3f}')")) =="
