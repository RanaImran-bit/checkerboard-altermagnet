#!/bin/bash
# run113.sh -- run_sync.sh adapted for the 113new SLURM cluster.
# Differences from the 250-series version:
#   * oneapi comes from `module load`, not /opt/intel/oneapi/setvars.sh
#   * python is the shared cluster anaconda (numpy 2.1.3 / scipy / pandas)
#   * everything stays under ~/Imran
# Usage: bash run113.sh <L> <U> <delta> <n> [NP] [t']
[ $# -ge 4 ] || { echo "usage: bash run113.sh <L> <U> <delta> <n> [NP] [t']"; exit 2; }
L=$1; U=$2; DELTA=$3; NFILL=$4; NP=${5:-8}; TPRIME=${6:-0.3}
T0=-1.0; V=0.0

source /etc/profile.d/modules.sh 2>/dev/null || true
module load oneapi/2022.3 2>/dev/null
set -u
PY=/opt/data/hpc/common/anaconda3/bin/python
SRC=$HOME/Imran/src
RUNS=$HOME/Imran/runs

T2=$($PY -c "print(f'{-1.0*$DELTA:.4f}')")
N=$($PY -c "print(round($NFILL*$L*$L/2))")
FOLD="$RUNS/L${L}n$(printf '%.3f' $NFILL)u${U}tA${T2}tt${TPRIME}N${N}"

[ $((L % 2)) -eq 0 ] || { echo "ERROR: L must be EVEN"; exit 1; }
command -v mpiifort >/dev/null 2>&1 || { echo "ERROR: mpiifort not found"; exit 1; }

if [ -f "$FOLD/nohup.out" ] && tail -3 "$FOLD/nohup.out" 2>/dev/null | grep -q time_sec; then
    echo "SKIP (done) $(basename $FOLD)"; exit 0
fi
mkdir -p "$FOLD"
for f in cp.f90 cpOut.f90 cpPara.f90 jiekou.f90 libuph.f90 mc2duph.f90 \
         parameter.f90 tk_check.f90 Makefile in.dat wf_unified.py; do
    cp "$SRC/$f" "$FOLD/"
done

sed -i.bak -E "s/parameter \(lx=[0-9]+,ly=[0-9]+/parameter (lx=${L},ly=${L}/" "$FOLD/parameter.f90"
sed -i.bak -E "s/parameter \(NUP=[0-9]+,NDN=[0-9]+/parameter (NUP=${N},NDN=${N}/" "$FOLD/parameter.f90"

$PY - "$FOLD/in.dat" <<PYEOF
import sys
fp=sys.argv[1]; lines=open(fp).readlines()
lines[4]=f"{$T0:.8f},{$TPRIME:.8f},{$T2:.8f},{0.0:.8f},-1\n"
lines[5]=f"{$U:.8f},{$V:.8f}\n"
open(fp,'w').writelines(lines)
PYEOF

( cd "$FOLD"
  $PY wf_unified.py > wf.log 2>&1 || { echo "TRIAL WF FAILED"; tail -3 wf.log; exit 1; }
  make clear >/dev/null 2>&1; make all > build.log 2>&1
  [ -f CPMC.exe ] || { echo "BUILD FAILED"; tail -5 build.log; exit 1; }
  echo "RUN $(basename $FOLD)  $(date +%H:%M)"
  mpirun -np $NP ./CPMC.exe > nohup.out 2>&1
  echo "DONE $(basename $FOLD)  $(date +%H:%M)"
)
