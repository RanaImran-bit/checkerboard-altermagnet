#!/bin/bash
# h-scan at L = 14 and L = 16, extending the L = 8/10/12 series.
#
# The h -> 0 intercept is nonzero at every size measured so far but falls with L
# (0.0327, 0.0236, 0.0184). The observed drop of 0.56x sits between the 1/L form
# (0.67, extrapolating to zero) and the 1/L^2 form (0.44, extrapolating to a finite
# value), so three points cannot fix the exponent. L = 16 is what discriminates:
# the two forms predict 0.0110 vs 0.0136 there, a 24% gap, against only 10% at
# L = 14. L = 14 fills the curve so the fit does not rest on the endpoints.
#
# Periodic boundaries and identical settings to the existing three sizes, so the
# five points form one series.
go () {
  h=$1; L=$2; nup=$3; np=$4
  ssh -n -o ConnectTimeout=20 -o BatchMode=yes -o ControlMaster=no "$h" \
    "cd ~/hscan && setsid nohup env L=$L NUPS=$nup US=4 DELTAS=0.3 \
     HS=0.05,0.1,0.2,0.3,0.5 NSEED=6 NPROC=$np \
     NW=320 NEQ=160 NBLK=30 BP=12 DT=0.05 \
     OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
     /opt/anaconda3/bin/python -u checkerboard_polarization.py > hscan_L${L}.log 2>&1 < /dev/null & sleep 2; echo '$h L=$L started'" &
}
go 250 14  98 30      # 196 sites, ~9 h
go 256 16 128 30      # 256 sites, ~20 h, the long pole
wait
