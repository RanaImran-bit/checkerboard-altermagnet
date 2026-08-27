#!/bin/bash
mkdir -p ~/Desktop/checkerboard-altermagnet/cb_model/checkerboard_work/data/mag_data
for N in 250 253 254 255 251; do
  scp "$N:~/mag_run/mag_L*_U*.csv"  ~/Desktop/checkerboard-altermagnet/cb_model/checkerboard_work/data/mag_data/ 2>/dev/null
  scp "$N:~/mag_run/sq_L*_U*.npz"   ~/Desktop/checkerboard-altermagnet/cb_model/checkerboard_work/data/mag_data/ 2>/dev/null
done
echo "collected to ~/Desktop/checkerboard-altermagnet/cb_model/checkerboard_work/data/mag_data/"
echo "csv: $(ls ~/Desktop/checkerboard-altermagnet/cb_model/checkerboard_work/data/mag_data/mag_L*_U*.csv 2>/dev/null | wc -l)  npz: $(ls ~/Desktop/checkerboard-altermagnet/cb_model/checkerboard_work/data/mag_data/sq_L*_U*.npz 2>/dev/null | wc -l)"
