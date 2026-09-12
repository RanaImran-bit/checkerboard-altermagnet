#!/bin/bash
# Node 256 -- first use. Two purposes:
#   1. smoke-test the node (conda path, detached launch, real throughput)
#   2. verify the fleet dashboard reports load correctly, by taking a KNOWN 8 cores
# The work is real: L=8 magnetic U=2, one of the gaps needed to complete the
# finite-size series. NPROC=8 leaves headroom per the group guide (256 has 11 free).
S=~/Desktop/checkerboard-altermagnet/code
ssh 256 'mkdir -p ~/mag_L8'
scp "$S/cpqmc.py" "$S/checkerboard.py" "$S/unified_scan.py" "$S/checkerboard_mag_scan.py" 256:~/mag_L8/
ssh -f 256 'cd ~/mag_L8 && SIZES=8 US=2 DELTAS=0,0.1,0.2,0.3,0.4 NPROC=8 NSEED=6 nohup /opt/anaconda3/bin/python checkerboard_mag_scan.py </dev/null > mag_L8.log 2>&1 &'
sleep 2; echo "== 256: magnetic L=8, U=2, delta 0-0.4, NPROC=8 =="
