#!/bin/bash
# node 251: PAIRING cube, L=12, U=6, high delta 0.5/0.6/0.7. All four channels.
# Separate dir ~/pair_hi/ so it cannot touch the delta<=0.4 pairing data.
# Driver writes after EACH delta, so a reboot costs one delta not the whole block.
S=~/Desktop/checkerboard-altermagnet/code
ssh 251 'mkdir -p ~/pair_hi'
scp "$S/cpqmc.py" "$S/checkerboard.py" "$S/checkerboard_fss_full.py" 251:~/pair_hi/
ssh -f 251 'cd ~/pair_hi && SIZES=12 US=6 DELTAS=0.5,0.6,0.7 NPROC=32 NSEED=6 nohup ~/anaconda3/bin/python checkerboard_fss_full.py </dev/null > pair_hi.log 2>&1 &'
sleep 1; echo "== 251: pairing L=12, U=6, delta 0.5/0.6/0.7 =="
