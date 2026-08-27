#!/bin/bash
# node 252 (32 physical cores, idle): dense-filling pairing scan, L=12, U=4
# PASS 1 = EVEN nup only (36,38,...,72) -> 19 fillings, Delta n = 0.028, ~16 h.
# The odd fillings are pass 2 on the main fleet once the magnetic run frees it.
S=~/Desktop/checkerboard-altermagnet/code
ssh 252 'mkdir -p ~/dense_run'
scp "$S/cpqmc.py" "$S/checkerboard.py" "$S/checkerboard_dense_n.py" 252:~/dense_run/
ssh -f 252 'cd ~/dense_run && L=12 U=4 NUP_MIN=36 NUP_MAX=72 STEP=2 NPROC=32 NSEED=6 nohup /opt/anaconda3/bin/python checkerboard_dense_n.py </dev/null > dense.log 2>&1 &'
sleep 1; echo "== 252 launched: dense pass 1, EVEN nup 36..72 (19 fillings, n=0.500..1.000) =="
