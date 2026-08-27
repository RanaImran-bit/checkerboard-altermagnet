#!/bin/bash
# PASS 2: the ODD nup (37,39,...,71) = 18 fillings, split over the main fleet.
# Run this ONLY after the magnetic run has finished on 250/253/254/255/251.
S=~/Desktop/checkerboard-altermagnet/code
i=0
for spec in "250:37:43:~/anaconda3/bin/python" "253:45:51:/opt/anaconda3/bin/python" \
            "254:53:59:/opt/anaconda3/bin/python" "255:61:65:/opt/anaconda3/bin/python" \
            "251:67:71:~/anaconda3/bin/python"; do
  N=${spec%%:*}; r=${spec#*:}; A=${r%%:*}; r=${r#*:}; B=${r%%:*}; PY=${r#*:}
  ssh $N 'mkdir -p ~/dense_run'
  scp "$S/cpqmc.py" "$S/checkerboard.py" "$S/checkerboard_dense_n.py" $N:~/dense_run/
  ssh -f $N "cd ~/dense_run && L=12 U=4 NUP_MIN=$A NUP_MAX=$B STEP=2 NPROC=32 NSEED=6 nohup $PY checkerboard_dense_n.py </dev/null > dense.log 2>&1 &"
  echo "== $N launched: dense pass 2, odd nup $A..$B =="
done
