#!/bin/bash
ssh -o BatchMode=yes 256 'echo "our python procs: $(ps -u phd25imran -o comm=|grep -c python)"; echo "load: $(cut -d" " -f1-3 /proc/loadavg)"; echo "rows: $(tail -n +2 ~/mag_L8/mag_L8_U2.csv 2>/dev/null|wc -l|tr -d " ")/180"; tail -2 ~/mag_L8/mag_L8.log 2>/dev/null'
