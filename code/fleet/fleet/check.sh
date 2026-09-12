#!/bin/bash
for N in 250 253 254 255 251; do
  echo "=== $N ==="
  ssh $N 'ls ~/cb_run/fss_full_L*_U*.csv 2>/dev/null | wc -l | xargs echo "blocks done:"; tail -2 ~/cb_run/fss_full.log'
done
