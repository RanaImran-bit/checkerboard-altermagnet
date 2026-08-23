#!/bin/bash
# Pairing VERTEX at L = 8 and L = 10, matched to the committed L = 12 run
# (beta = 32, N_w = 500, tau window 0.8, 6 seeds, all seven delta).
#
# Closes the one real gap in the paper's strongest claim: the vertex, and with it
# the d_xy dominance for delta >= 0.3, currently exists at L = 12 only. The
# multi-size data already on disk is FULL chi, which is dominated by the
# uncorrelated bubble and therefore cannot settle the question.
#
# U = 4 is representative: it sits inside the range where d_xy leads at L = 12
# and above the delta ~ 0.3 channel crossover.
go () {
  h=$1; L=$2; nup=$3; np=$4
  ssh -n -o ConnectTimeout=20 -o BatchMode=yes -o ControlMaster=no "$h" \
    "cd ~/chi_fss && setsid nohup env L=$L NUPS=$nup US=4 DELTAS=0.1,0.2,0.3,0.4,0.5,0.6,0.7 \
     NSEED=6 NPROC=$np NW=500 NEQ=640 NBLK=40 BP=16 DT=0.05 \
     OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
     /opt/anaconda3/bin/python -u checkerboard_eqtime.py > chifss_L${L}.log 2>&1 < /dev/null & sleep 2; echo '$h L=$L started'" &
}
go 252 8  32 42     # 64 sites,  ~1 h per cell
go 251 10 50 42     # 100 sites, ~4 h per cell, sets the wall time
wait
