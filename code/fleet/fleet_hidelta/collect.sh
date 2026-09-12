#!/bin/bash
mkdir -p ~/Desktop/checkerboard-altermagnet/cb_model/checkerboard_work/data/mag_hi_data
for N in 255 252 254 250 251; do
  scp -q "$N:~/mag_hi/mag_L12_U*.csv" ~/Desktop/checkerboard-altermagnet/cb_model/checkerboard_work/data/mag_hi_data/ 2>/dev/null
  scp -q "$N:~/mag_hi/sq_L12_U*.npz"  ~/Desktop/checkerboard-altermagnet/cb_model/checkerboard_work/data/mag_hi_data/ 2>/dev/null
done
echo "collected: $(ls ~/Desktop/checkerboard-altermagnet/cb_model/checkerboard_work/data/mag_hi_data/mag_L12_U*.csv 2>/dev/null | wc -l) of 5"
