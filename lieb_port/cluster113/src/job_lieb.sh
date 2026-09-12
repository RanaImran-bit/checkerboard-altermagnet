#!/bin/bash
#SBATCH -J lieb_am
#SBATCH -p all
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH -t 30-00:00:00
#SBATCH --mem=16G
#SBATCH -o slurm-%j.out
# 113new has no Intel oneAPI for this account, so the GNU Fortran + OpenMPI
# Makefile is used and the trial wave function needs anaconda's python: the
# default python3 on this cluster has no numpy.
# Intel oneAPI, NOT the GNU/OpenMPI stack. Measured on this cluster: the GNU
# build against reference LAPACK/BLAS ran one measurement block per 15 minutes
# on 12 sites, about 60x slower than the same job on 251. ifort + MKL is what
# 250-255 use and it is available here too.
# One thread per MPI rank. Without this each rank spawns a full MKL/OpenMP
# thread pool sized to the node, and n001 reached a load of 652 on 96 cores with
# only 80 ranks running. That oversubscription, not the arithmetic, was what
# made the first round of jobs look impossibly slow.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
source /opt/data/hpc/common/intel/oneapi/setvars.sh >/dev/null 2>&1
cp Makefile.intel Makefile
$HOME/anaconda3/bin/python3 wf_lieb.py > wf.log 2>&1 \
  || { echo "TRIAL WF FAILED"; tail -20 wf.log; exit 1; }
make clear >/dev/null 2>&1
make all > build.log 2>&1
[ -f CPMC.exe ] || { echo "BUILD FAILED"; tail -25 build.log; exit 1; }
mpirun -np ${SLURM_NTASKS} ./CPMC.exe
