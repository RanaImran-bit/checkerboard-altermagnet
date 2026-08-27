#!/bin/bash
# launches the whole magnetic fleet (5 nodes, split by U). Expect ~5 password prompts x3.
for s in launch_mag_250.sh launch_mag_253.sh launch_mag_254.sh launch_mag_255.sh launch_mag_251.sh; do
  bash ~/Desktop/fleet_mag/$s
done
