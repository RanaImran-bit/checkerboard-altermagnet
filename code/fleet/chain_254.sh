#!/bin/bash
# Wait for the delta<=0.4 job on this node to exit, then run the delta=0.5,0.6,0.7 work.
# Output goes to FRESH dirs: the drivers name files mag_L{L}_U{U}.csv with no delta in
# the name, so writing into the old dirs would silently overwrite the delta<=0.4 data.
while pgrep -f 'checkerboard_mag_scan\.py' > /dev/null; do sleep 60; done
echo "$(date '+%F %T') predecessor finished, starting magnetic" >> $HOME/chain_254.log

mkdir -p $HOME/mag_hi2 && cd $HOME/mag_hi2
cp $HOME/mag_gap1/checkerboard_mag_scan.py $HOME/mag_gap1/checkerboard.py \
   $HOME/mag_gap1/cpqmc.py $HOME/mag_gap1/unified_scan.py .
SIZES=10 US=6,8 DELTAS=0.5,0.6,0.7 NPROC=32 NSEED=6 \
  /opt/anaconda3/bin/python checkerboard_mag_scan.py < /dev/null > mag_hi2.log 2>&1
echo "$(date '+%F %T') magnetic done, starting pairing" >> $HOME/chain_254.log

cd $HOME/pair_hi2
SIZES=10 US=6,8 DELTAS=0.5,0.6,0.7 NPROC=32 NSEED=6 \
  /opt/anaconda3/bin/python checkerboard_fss_full.py < /dev/null > pair_hi2.log 2>&1
echo "$(date '+%F %T') pairing done, node idle" >> $HOME/chain_254.log
