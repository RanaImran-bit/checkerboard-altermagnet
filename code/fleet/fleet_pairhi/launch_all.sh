#!/bin/bash
# Run ONLY after the magnetic high-delta run has finished (check with fleet_hidelta/check.sh)
for s in launch_pairhi_255.sh launch_pairhi_252.sh launch_pairhi_254.sh launch_pairhi_250.sh launch_pairhi_251.sh; do
  bash ~/Desktop/checkerboard-altermagnet/cb_model/checkerboard_work/fleet/fleet_pairhi/$s
done
