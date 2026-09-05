#!/bin/bash
# =============================================================================
# campaign.sh -- L=14 checkerboard campaign with every known fix applied.
#
#   bash campaign.sh build   <L> <N>                 compile one binary
#   bash campaign.sh trial   <L> <U> <delta> <N>     generate one cell's trial
#   bash campaign.sh run     <cell> <seedidx>        launch one 1-core run
#   bash campaign.sh cells   <L> <N>                 list the 49 cells of a block
#
# Fixes carried by this campaign, all validated earlier:
#   * trial chosen by the VARIATIONAL energy, not the eigenvalue sum. The old
#     criterion picked a worse trial in 36% of 304 audited runs and, at half
#     filling, always overstated the staggered moment.
#   * SCF tightened to 1e-8 with 4000 iterations and 12 seeds, because fresh
#     solves were beating the saved trials.
#   * unique per-run seed via seed.dat. The clock seed spans only milliseconds,
#     so one-core runs launched together would otherwise share seeds and report
#     error bars that are too small.
#   * progress flushed, so a running job is visible instead of silent for hours.
#
# One core per run: the binary is launched with mpirun -np 1, so a 64-core box
# carries 64 independent runs. The trial depends on the cell but NOT on the seed,
# so it is generated once per cell and shared by that cell's seeds. parameter.f90
# compiles in L and N only, so one binary serves all U and delta of a block.
# =============================================================================
set -u
ROOT="$HOME/campaign_L14L16"
SRC="$ROOT/src"; BIN="$ROOT/bin"; TRI="$ROOT/trials"; RUN="$ROOT/runs"
TPRIME=0.3; T0=-1.0; V=0.0
PY=/opt/anaconda3/bin/python3; [ -x "$PY" ] || PY="$HOME/anaconda3/bin/python3"
setvars(){ [ -f /opt/intel/oneapi/setvars.sh ] && { set +u; . /opt/intel/oneapi/setvars.sh >/dev/null 2>&1; set -u; }; }
cellname(){ # L U delta N
  local t2; t2=$($PY -c "print(f'{-1.0*$3:.4f}')")
  local nf; nf=$($PY -c "print(f'{2*$4/($1*$1):.3f}')")
  echo "L$1n${nf}u$2tA${t2}tt${TPRIME}N$4"; }

case "${1:-}" in
build)
  L=$2; N=$3; setvars; mkdir -p "$BIN" "$BIN/.w_${L}_${N}"
  cd "$BIN/.w_${L}_${N}" && cp "$SRC"/*.f90 "$SRC"/Makefile .
  sed -i -E "s/parameter \(lx=[0-9]+,ly=[0-9]+/parameter (lx=${L},ly=${L}/" parameter.f90
  sed -i -E "s/parameter \(NUP=[0-9]+,NDN=[0-9]+/parameter (NUP=${N},NDN=${N}/" parameter.f90
  make clear >/dev/null 2>&1; make all > build.log 2>&1
  [ -f CPMC.exe ] || { echo "BUILD FAILED L=$L N=$N"; grep -i error build.log | head -5; exit 1; }
  mv CPMC.exe "$BIN/CPMC_L${L}_N${N}.exe"; echo "built bin/CPMC_L${L}_N${N}.exe" ;;

trial)
  L=$2; U=$3; DE=$4; N=$5; C=$(cellname "$L" "$U" "$DE" "$N"); D="$TRI/$C"
  [ -f "$D/wfup.txt" ] && { echo "trial exists: $C"; exit 0; }
  mkdir -p "$D"; cp "$SRC/in.dat" "$SRC/wf_unified.py" "$D/"
  T2=$($PY -c "print(f'{-1.0*$DE:.4f}')")
  $PY - "$D/in.dat" <<PYEOF
import sys
fp=sys.argv[1]; L=open(fp).readlines()
L[4]=f"{$T0:.8f},{$TPRIME:.8f},{$T2:.8f},{0.0:.8f},-1\n"
L[5]=f"{$U:.8f},{$V:.8f}\n"
open(fp,'w').writelines(L)
PYEOF
  ( cd "$D" && OMP_NUM_THREADS=1 $PY wf_unified.py > wf.log 2>&1 ) \
    && echo "trial ok: $C" || { echo "TRIAL FAILED: $C"; tail -3 "$D/wf.log"; exit 1; } ;;

run)
  C=$2; K=$3; L=$(echo "$C" | sed -E "s/^L([0-9]+)n.*/\1/"); N=$(echo "$C" | sed -E "s/.*N([0-9]+)$/\1/")
  D="$RUN/${C}_s${K}"
  [ -f "$D/corr.dat" ] && { echo "done already: ${C}_s${K}"; exit 0; }
  mkdir -p "$D"; cp "$TRI/$C/in.dat" "$TRI/$C/wfup.txt" "$TRI/$C/wfdn.txt" "$D/" || exit 1
  cp "$BIN/CPMC_L${L}_N${N}.exe" "$D/CPMC.exe" || exit 1
  # unique seed: cell hash mixed with the seed index, so no two runs collide
  $PY -c "import zlib;print(abs(zlib.crc32(b'$C'))%900000000 + $K*7919 + 1)" > "$D/seed.dat"
  setvars
  ( cd "$D" && nohup mpirun -np 1 ./CPMC.exe > nohup.out 2>&1 & )
  echo "launched ${C}_s${K}" ;;

runfg)  # same as run but blocks, so xargs -P can cap concurrency at the core count
  C=$2; K=$3; L=$(echo "$C" | sed -E "s/^L([0-9]+)n.*/\1/"); N=$(echo "$C" | sed -E "s/.*N([0-9]+)$/\1/")
  D="$RUN/${C}_s${K}"
  [ -f "$D/corr.dat" ] && { echo "skip done ${C}_s${K}"; exit 0; }
  mkdir -p "$D"; cp "$TRI/$C/in.dat" "$TRI/$C/wfup.txt" "$TRI/$C/wfdn.txt" "$D/" || exit 1
  cp "$BIN/CPMC_L${L}_N${N}.exe" "$D/CPMC.exe" || exit 1
  $PY -c "import zlib;print(abs(zlib.crc32(b'$C'))%900000000 + $K*7919 + 1)" > "$D/seed.dat"
  setvars
  ( cd "$D" && OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 mpirun -np 1 ./CPMC.exe > nohup.out 2>&1 )
  echo "finished ${C}_s${K}" ;;

cells)
  L=$2; N=$3
  for U in 0.0 2.0 3.0 3.5 4.0 4.5 5.0; do for DE in 0.1 0.2 0.3 0.4 0.5 0.6 0.7; do
    echo "$(cellname "$L" "$U" "$DE" "$N") $U $DE"; done; done ;;
*) echo "usage: campaign.sh {build|trial|run|runfg|cells} ..."; exit 2 ;;
esac
