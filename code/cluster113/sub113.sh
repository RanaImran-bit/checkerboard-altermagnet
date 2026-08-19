#!/bin/bash
#SBATCH -p all
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -o %x-%j.out
# Partition MUST be `all`: the default `normal` partition has MaxTime=00:01:00
# and would kill every job after 60 seconds.
bash $HOME/Imran/run113.sh "$@"
