#!/bin/bash
# Submit ONE Lieb-lattice CPQMC job to the 113new Slurm cluster.
#
#   bash submit_lieb.sh <cells> <U> <t2> <n_cell>
#
# Same physics and same in.dat slot mapping as run_lieb.sh on 250-255; the only
# differences are that the build and the run happen inside the batch job, the
# GNU/OpenMPI Makefile is used, and srun replaces mpirun.
#
# Partition MUST be "all": "normal" has a one hour limit and these run for days.
set -eu
[ $# -ge 4 ] || { echo "usage: bash submit_lieb.sh <cells> <U> <t2> <n_cell>"; exit 2; }
CELLS=$1; U=$2; T2=$3; NCELL=$4
EPS_A=${EPS_A:-0.0}; T0=-1.0; V=0.0
PY=$HOME/anaconda3/bin/python3
SRC="$(cd "$(dirname "$0")/src" && pwd)"

L=$((2*CELLS)); NSITE=$((3*CELLS*CELLS)); NELEC=$((NCELL*CELLS*CELLS))
if [ "$NCELL" -eq 3 ]; then                 # half filling: Lieb's theorem, ferrimagnetic
  NUP=$((2*CELLS*CELLS)); NDN=$((CELLS*CELLS))
else                                        # compensated altermagnet
  [ $((NELEC % 2)) -eq 0 ] || { echo "ERROR: NELEC=$NELEC odd"; exit 1; }
  NUP=$((NELEC/2)); NDN=$NUP
fi
ET=${ETRIAL:-$($PY -c "print(f'{-1.0*$NSITE:.6f}')")}
NFILL=$($PY -c "print(f'{$NELEC/$NSITE:.3f}')")
FOLD="L${L}n${NFILL}u${U}tA${T2}tt${T0}N${NUP}"

[ -d "$FOLD" ] && { echo "SKIP (exists): $FOLD"; exit 0; }
mkdir "$FOLD"; cp "$SRC"/* "$FOLD/"

sed -i -E "s/parameter \(lx=[0-9]+,ly=[0-9]+/parameter (lx=${L},ly=${L}/" "$FOLD/parameter.f90"
sed -i -E "s/parameter \(NUP=[0-9]+,NDN=[0-9]+/parameter (NUP=${NUP},NDN=${NDN}/" "$FOLD/parameter.f90"
$PY - "$FOLD/in.dat" <<PYEOF
import sys
fp=sys.argv[1]; L=open(fp).readlines()
L[2]=L[2].split(",")[0] + f",{$ET:.6f}\n"
L[4]=f"{$T0:.8f},{$T2:.8f},{$EPS_A:.8f},{0.0:.8f},-1\n"
L[5]=f"{$U:.8f},{$V:.8f}\n"
open(fp,'w').writelines(L)
PYEOF

# NTASKS overrides the rank count baked into job_lieb.sh. Small systems can fail
# with too many ranks ("Problem in sgefa routine"), so it is tunable per run.
NTASKS=${NTASKS:-16}
JID=$(cd "$FOLD" && sbatch --parsable --ntasks-per-node=$NTASKS job_lieb.sh)
echo "submitted $JID  $FOLD   sites=$NSITE NUP=$NUP NDN=$NDN ranks=$NTASKS U(B,C)=$U t2=$T2 eps_A=$EPS_A"
