#!/bin/bash
for N in 250 253 254 255 251; do
  echo "=== $N ==="
  ssh $N 'ls ~/mag_run/mag_L*_U*.csv 2>/dev/null | wc -l | xargs echo "blocks done (of 3):"; tail -2 ~/mag_run/mag.log'
done
