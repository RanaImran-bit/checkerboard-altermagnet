#!/bin/bash
# pull all finished per-(L,U) CSVs from every node into ~/Desktop/checkerboard-altermagnet/cb_model/checkerboard_work/data/fss_data/
mkdir -p ~/Desktop/checkerboard-altermagnet/cb_model/checkerboard_work/data/fss_data
for N in 250 253 254 255 251; do scp "$N:~/cb_run/fss_full_L*_U*.csv" ~/Desktop/checkerboard-altermagnet/cb_model/checkerboard_work/data/fss_data/ 2>/dev/null; done
echo "collected to ~/Desktop/checkerboard-altermagnet/cb_model/checkerboard_work/data/fss_data/"; ls ~/Desktop/checkerboard-altermagnet/cb_model/checkerboard_work/data/fss_data/fss_full_L*_U*.csv 2>/dev/null | wc -l | xargs echo "per-(L,U) files:"
