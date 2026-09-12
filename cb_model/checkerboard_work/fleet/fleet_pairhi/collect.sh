#!/bin/bash
D=~/Desktop/checkerboard-altermagnet/cb_model/checkerboard_work/data/pair_hi_data
mkdir -p $D
for N in 255 250 251 253; do scp -q "$N:~/pair_hi/fss_full_L12_U*.csv" $D/ 2>/dev/null; done
echo "collected: $(ls $D/fss_full_L12_U*.csv 2>/dev/null | wc -l) of 5 U blocks"
