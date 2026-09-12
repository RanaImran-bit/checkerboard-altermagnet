#!/bin/bash
while pgrep -f "checkerboard_fss_full\.py" > /dev/null; do sleep 60; done
echo "$(date "+%F %T") delta-extension finished, starting pairfull L=10" >> $HOME/chain_254.log
cd $HOME/pairfull
SIZES=10 US=0,2,4,6,8 DELTAS=0 NPROC=32 NSEED=6   /opt/anaconda3/bin/python checkerboard_pairfull.py < /dev/null > pairfull.log 2>&1
echo "$(date "+%F %T") pairfull L=10 done" >> $HOME/chain_254.log
