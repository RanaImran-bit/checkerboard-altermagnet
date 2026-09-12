#!/bin/bash
mkdir -p ~/Desktop/checkerboard-altermagnet/cb_model/checkerboard_work/data/pair_hi_data
for N in 255 252 254 250 251; do scp -q "$N:~/pair_hi/fss_full_L12_U*.csv" ~/Desktop/checkerboard-altermagnet/cb_model/checkerboard_work/data/pair_hi_data/ 2>/dev/null; done
echo "collected: $(ls ~/Desktop/checkerboard-altermagnet/cb_model/checkerboard_work/data/pair_hi_data/fss_full_L12_U*.csv 2>/dev/null | wc -l) of 5"
