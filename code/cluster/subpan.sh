#!/bin/bash
#SBATCH -p all
#SBATCH -N 1
#SBATCH -n 18
#SBATCH --mem-per-cpu=12G
#SBATCH -o logs/%x-%j.out
# Pan et al. channel comparison at the ONE point where the two Hamiltonians
# coincide: delta = 0.3, half filling. Their basis (s, d) against ours
# (on-site s, ext s, d_x2-y2, d_xy), same walkers, same beta, same seeds.
B=/opt/data/hpc/ZONES/all/PROJECTS/a1/PRIVATE/jianyu/Imran/checkerboard/pan
cd "$B" || exit 1
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export L=12 NUPS=72 DELTAS=0.3 US=3,4,4.5 NSEED=6 NPROC=18
export NW=500 NEQ=640 NBLK=40 BP=16 DT=0.05
/opt/data/hpc/common/anaconda3/bin/python -u pan_compare.py
