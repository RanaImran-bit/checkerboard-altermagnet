#!/bin/bash
# One U per node, 6 seeds each: 6 concurrent cells per node instead of 18, so
# memory contention stays low and each cell runs near the 11.3 h single-process
# time measured on this hardware. Wall time = ONE cell, since parallelism across
# cells is already saturated at one cell per core -- more nodes cannot shorten it
# further, only faster cores and lower contention can.
go () {
  h=$1; U=$2
  ssh -n -o ConnectTimeout=20 -o BatchMode=yes -o ControlMaster=no "$h" \
    "cd ~/pan_compare && setsid nohup env L=12 NUPS=72 DELTAS=0.3 US=$U NSEED=6 NPROC=6 \
     NW=500 NEQ=640 NBLK=40 BP=16 DT=0.05 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
     /opt/anaconda3/bin/python -u pan_compare.py > pan_U${U}.log 2>&1 < /dev/null & sleep 2; echo '$h U=$U started'" &
}
go 252 3
go 253 4
go 256 4.5
wait
