#!/bin/bash
for N in 252 254; do echo -n "$N: "; ssh -o BatchMode=yes -o ConnectTimeout=8 $N 'echo -n "procs=$(ps -u phd25imran -o comm=|grep -c python) rows=$(tail -n +2 ~/mag_hi/mag_L12_U*.csv 2>/dev/null|wc -l|tr -d " ") | "; tail -1 ~/mag_hi/mag_hi.log 2>/dev/null|cut -c1-55' 2>&1|head -1; done
