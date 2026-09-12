#!/bin/bash
mkdir -p ~/Desktop/checkerboard-altermagnet/cb_model/checkerboard_work/data/dense_data
for N in 252 250 253 254 255 251; do scp "$N:~/dense_run/dense_L12_U4_nup*.csv" ~/Desktop/checkerboard-altermagnet/cb_model/checkerboard_work/data/dense_data/ 2>/dev/null; done
for N in 255 251; do scp "$N:~/mag_run/mag_L12_U*.csv" ~/Desktop/checkerboard-altermagnet/cb_model/checkerboard_work/data/dense_data/ 2>/dev/null
                    scp "$N:~/mag_run/sq_L12_U*.npz"  ~/Desktop/checkerboard-altermagnet/cb_model/checkerboard_work/data/dense_data/ 2>/dev/null; done
echo "dense: $(ls ~/Desktop/checkerboard-altermagnet/cb_model/checkerboard_work/data/dense_data/dense_*.csv 2>/dev/null|wc -l)  mag: $(ls ~/Desktop/checkerboard-altermagnet/cb_model/checkerboard_work/data/dense_data/mag_*.csv 2>/dev/null|wc -l)"
