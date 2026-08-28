#!/bin/bash
# =============================================================================
# run_one.sh — launch ONE checkerboard CPQMC job on an Intel-turnkey machine
#              (thkclusters ports 250,251,252,253,254,255 — oneAPI + MKL).
#
#   Run from INSIDE ~/Checkerboard_Model/ :
#     bash run_one.sh <L> <U> <delta> <n> [NP]
#
#   examples:
#     bash run_one.sh 18 4.0 0.3 1.0 8      # finite-size L=18, U=4, δ=0.3, half filling
#     bash run_one.sh 16 4.0 0.3 1.0 8      # finite-size L=16, δ=0.3
#     bash run_one.sh 14 4.0 0.2 0.847 8    # Phase-3 doping (n=0.847), etc.
#
#   It: sets the model params, builds with the Intel Makefile, and launches
#   mpirun -np NP ./CPMC.exe  DETACHED (nohup) so it survives logout.
#   Fixed model convention (paper 2605.11669): t=1 -> t0=-1 ; t'=-0.3 -> t1=+0.3 ;
#   V=0 ; anisotropy delta -> in.dat t2 = -delta ; tam=0 (spin-independent).
#   L MUST be even. n=1 is half filling (N=L^2/2).
# =============================================================================
set -u
[ $# -ge 4 ] || { echo "usage: bash run_one.sh <L> <U> <delta> <n> [NP]"; exit 2; }
L=$1; U=$2; DELTA=$3; NFILL=$4; NP=${5:-8}

T0=-1.0; TPRIME=0.3; V=0.0                 # fixed model params

# Intel toolchain (mpiifort + MKL)
[ -f /opt/intel/oneapi/setvars.sh ] && { set +u; source /opt/intel/oneapi/setvars.sh >/dev/null 2>&1; set -u; }

# python with numpy (for the trial WF)
PY=python3
for p in /opt/anaconda3/bin/python3 "$HOME/anaconda3/bin/python3"; do [ -x "$p" ] && PY="$p"; done

T2=$($PY -c "print(f'{-1.0*$DELTA:.4f}')")          # t2 = -delta (raw element)
N=$($PY -c "print(round($NFILL*$L*$L/2))")          # N_up=N_dn = round(n*L^2/2)
FOLD="L${L}n$(printf '%.3f' $NFILL)u${U}tA${T2}tt${TPRIME}N${N}"

echo ">>> $FOLD   host=$(hostname)  NP=$NP  py=$PY"
[ $((L % 2)) -eq 0 ] || { echo "ERROR: L must be EVEN for the checkerboard"; exit 1; }
command -v mpiifort >/dev/null 2>&1 || { echo "ERROR: mpiifort not found (source /opt/intel/oneapi/setvars.sh?)"; exit 1; }

mkdir -p "$FOLD"
for f in *.f90 Makefile in.dat wf_unified.py; do cp "$f" "$FOLD/"; done

sed -i.bak -E "s/parameter \(lx=[0-9]+,ly=[0-9]+/parameter (lx=${L},ly=${L}/" "$FOLD/parameter.f90"
sed -i.bak -E "s/parameter \(NUP=[0-9]+,NDN=[0-9]+/parameter (NUP=${N},NDN=${N}/" "$FOLD/parameter.f90"

$PY - "$FOLD/in.dat" <<PYEOF
import sys
fp=sys.argv[1]; lines=open(fp).readlines()
lines[4]=f"{$T0:.8f},{$TPRIME:.8f},{$T2:.8f},{0.0:.8f},-1\n"   # t0,t1,t2,tam,alphat1
lines[5]=f"{$U:.8f},{$V:.8f}\n"                                # ud,vpd
open(fp,'w').writelines(lines)
PYEOF

( cd "$FOLD"
  $PY wf_unified.py > wf.log 2>&1 || { echo "TRIAL WF FAILED (see $FOLD/wf.log)"; exit 1; }
  make clear >/dev/null 2>&1; make all > build.log 2>&1
  [ -f CPMC.exe ] || { echo "BUILD FAILED (see $FOLD/build.log)"; tail -6 build.log; exit 1; }
  nohup mpirun -np $NP ./CPMC.exe > nohup.out 2>&1 &
  echo "STARTED pid=$!  ->  $FOLD/nohup.out   (done when last line is 'time_sec ...')"
)
