#!/bin/bash
for N in 252 254 255; do echo -n "$N: "; ssh -o BatchMode=yes -o ConnectTimeout=8 $N 'echo -n "procs=$(ps -u phd25imran -o comm=|grep -c python) rows=$(tail -n +2 ~/mag_gap1/mag_L*_U*.csv 2>/dev/null|wc -l|tr -d " ")/180 | "; grep "delta=" ~/mag_gap1/mag_gap1.log 2>/dev/null|tail -1|sed "s/^ *//"' 2>&1|head -1; done
