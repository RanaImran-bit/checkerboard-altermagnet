#!/bin/bash
D=~/Desktop/checkerboard-altermagnet/cb_model/checkerboard_work/data/mag_gap1_data
mkdir -p $D
for N in 252 254 255; do scp -q "$N:~/mag_gap1/mag_L*_U*.csv" $D/ 2>/dev/null; done
echo "collected: $(ls $D/*.csv 2>/dev/null|wc -l) of 3"
