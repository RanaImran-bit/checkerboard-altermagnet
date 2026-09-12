#!/bin/bash
# Launch the Lieb structure-factor runs, one node per (size, U) slice, all in
# parallel so no single slow ssh stalls the rest.
go () {
  h=$1; tag=$2; shift 2
  ssh -n -o ConnectTimeout=20 -o BatchMode=yes -o ControlMaster=no "$h" \
    "cd ~/lieb-altermagnet/code && setsid nohup env $* OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /opt/anaconda3/bin/python -u lieb_sq.py > sq_${tag}.log 2>&1 < /dev/null & sleep 2; echo '$h:$tag started'" &
}
go 250 L8u2-4-5   SIZES=8 US=2,4,5 NSEED=6 NPROC=18
go 252 L8u6-7     SIZES=8 US=6,7   NSEED=6 NPROC=12
go 254 L8u8       SIZES=8 US=8     NSEED=6 NPROC=6
go 255 L8u10      SIZES=8 US=10    NSEED=6 NPROC=6
wait
