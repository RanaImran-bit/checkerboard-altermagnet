#!/bin/bash
for N in 255 252 254 250 251; do
  echo -n "$N: "
  ssh -o BatchMode=yes -o ConnectTimeout=8 $N 'echo -n "procs=$(ps -u phd25imran -o comm=|grep -c python) csv=$(ls ~/pair_hi/fss_full_L12_U*.csv 2>/dev/null|wc -l|tr -d " ") | "; tail -1 ~/pair_hi/pair_hi.log 2>/dev/null | cut -c1-60' 2>&1 | head -1
done
