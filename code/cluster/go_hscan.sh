#!/bin/bash
# Delta_tot vs pinning field h at three sizes, for the h -> 0 extrapolation.
#
# THE decisive test for the manuscript. A finite intercept that survives growing L
# is spontaneous order; an intercept consistent with zero means the splitting is
# entirely field-induced and the paper must be reframed. S(pi,pi)/N already falls
# as 1/N at every delta and U up to 8, so this is the claim that needs settling.
#
# One size per node so the wall time is a single L=12 cell, not a queue. beta is
# raised from the exploratory default of 3 to 8 (NEQ=160) and walkers from 160 to
# 320, since an order-parameter intercept needs a longer projection than a
# boundary-condition comparison did.
go () {
  h=$1; L=$2; nup=$3; np=$4
  ssh -n -o ConnectTimeout=20 -o BatchMode=yes -o ControlMaster=no "$h" \
    "cd ~/hscan && setsid nohup env L=$L NUPS=$nup US=4 DELTAS=0.3 \
     HS=0.05,0.1,0.2,0.3,0.5 NSEED=6 NPROC=$np \
     NW=320 NEQ=160 NBLK=30 BP=12 DT=0.05 \
     OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
     /opt/anaconda3/bin/python -u checkerboard_polarization.py > hscan_L${L}.log 2>&1 < /dev/null & sleep 2; echo '$h L=$L started'" &
}
go 253 8  32 30      # 64 sites,  cheapest
go 256 10 50 30      # 100 sites
go 255 12 72 30      # 144 sites, sets the wall time
wait
