#!/bin/bash
# run_sync.sh -- run_one.sh but SYNCHRONOUS, so a worker loop can queue jobs.
# The original ends with `nohup mpirun ... &` which returns immediately; that is
# right for one-off launches but means a queue would start everything at once.
# Usage: bash run_sync.sh <L> <U> <delta> <n> [NP] [t']
[ $# -ge 4 ] || { echo "usage: bash run_sync.sh <L> <U> <delta> <n> [NP] [t']"; exit 2; }
L=$1; U=$2; DELTA=$3; NFILL=$4; NP=${5:-8}; TPRIME=${6:-0.3}
T0=-1.0; V=0.0

# oneapi's setvars.sh references unset variables, which is FATAL under set -u
# and kills the shell silently with status 0. Source it first, then turn on -u.
. /opt/intel/oneapi/setvars.sh >/dev/null 2>&1
set -u
PY=python3
for p in /opt/anaconda3/bin/python3 "$HOME/anaconda3/bin/python3"; do [ -x "$p" ] && PY="$p"; done

T2=$($PY -c "print(f'{-1.0*$DELTA:.4f}')")
N=$($PY -c "print(round($NFILL*$L*$L/2))")
FOLD="L${L}n$(printf '%.3f' $NFILL)u${U}tA${T2}tt${TPRIME}N${N}"

[ $((L % 2)) -eq 0 ] || { echo "ERROR: L must be EVEN"; exit 1; }
command -v mpiifort >/dev/null 2>&1 || { echo "ERROR: mpiifort not found"; exit 1; }

# already finished? nohup.out ends with the timing line when complete
if [ -f "$FOLD/nohup.out" ] && tail -3 "$FOLD/nohup.out" 2>/dev/null | grep -q time_sec; then
    echo "SKIP (done) $FOLD"; exit 0
fi
mkdir -p "$FOLD"
for f in *.f90 Makefile in.dat wf_unified.py; do cp "$f" "$FOLD/"; done

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
  $PY wf_unified.py > wf.log 2>&1 || { echo "TRIAL WF FAILED $FOLD"; exit 1; }
  make clear >/dev/null 2>&1; make all > build.log 2>&1
  [ -f CPMC.exe ] || { echo "BUILD FAILED $FOLD"; tail -4 build.log; exit 1; }
  echo "RUN $FOLD  $(date +%H:%M)"
  mpirun -np $NP ./CPMC.exe > nohup.out 2>&1          # <-- FOREGROUND
  echo "DONE $FOLD  $(date +%H:%M)"
)
