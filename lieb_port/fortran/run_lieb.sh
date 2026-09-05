#!/bin/bash
# =============================================================================
# run_lieb.sh -- launch ONE Lieb-lattice CPQMC job.
#
#   bash run_lieb.sh <cells> <U> <t2> <n_cell> [NP]
#
#     cells  : unit cells per direction. The FINE grid is 2*cells, so the Lieb
#              lattice has 3*cells^2 sites and only those sites exist: the
#              remaining quarter of the fine grid is never indexed.
#     U      : Hubbard U, applied to sublattices B and C ONLY (A is the oxygen).
#     t2     : next-nearest B-C hopping. REQUIRED for the altermagnet; 0 keeps
#              the flat band. Goes into the in.dat t1 slot.
#     n_cell : ELECTRONS PER UNIT CELL, out of 6. This is the paper's n.
#                2, 4 -> altermagnetic Mott insulator, COMPENSATED, NUP = NDN
#                3    -> half filling, ferrimagnetic by Lieb's theorem
#
#   env: EPS_A on-site energy on sublattice A (default 0), ETRIAL, SEED_ONLY
#
#   examples:
#     bash run_lieb.sh 4 10.0 0.3 4 8     # 48 sites, the published AM point
#     bash run_lieb.sh 4 4.0  0.0 3 8     # half filling, the ferrimagnet
#
# Model: Kaushal and Franz, PRL, arXiv:2412.16421, Eq. (1). in.dat slots are
# REPURPOSED and the names do not match the physics:
#     in.dat t0  = nearest neighbour A-B and A-C   (the paper's t1 = -1)
#     in.dat t1  = next nearest B-C                (the paper's t2)
#     in.dat t2  = eps_A, on-site energy on A
#     in.dat tam = unused, must be 0
#
# NSTATES is 3*cells^2, the real sites and nothing else. The k mesh is still the
# fine grid's lx*ly points, four times the Lieb Brillouin zone, so k-space
# quantities are indexed by NKPTS and not by NSTATES.
# =============================================================================
set -u
[ $# -ge 4 ] || { echo "usage: bash run_lieb.sh <cells> <U> <t2> <n_cell> [NP]"; exit 2; }
CELLS=$1; U=$2; T2=$3; NCELL=$4; NP=${5:-8}
EPS_A=${EPS_A:-0.0}

T0=-1.0; V=0.0                       # NN amplitude, fixed at the paper's t1
L=$((2*CELLS)); NSITE=$((3*CELLS*CELLS)); NELEC=$((NCELL*CELLS*CELLS))

[ -f /opt/intel/oneapi/setvars.sh ] && { set +u; source /opt/intel/oneapi/setvars.sh >/dev/null 2>&1; set -u; }
PY=python3
for p in /opt/anaconda3/bin/python3 "$HOME/anaconda3/bin/python3"; do [ -x "$p" ] && PY="$p"; done

# The altermagnet is COMPENSATED: the moments sit on B and C, antialigned, and
# the total is zero. It therefore lives at NUP = NDN and CANNOT be reached from
# a sector with a net moment imposed by hand.
#
# Half filling is the exception. There Lieb's theorem applies (bipartite, with
# |A| = cells^2 and |BuC| = 2*cells^2) and forces S = cells^2/2, a ferrimagnet.
# Forcing NUP = NDN there also puts the Fermi level inside the cells^2-fold
# degenerate flat band, where the trial determinant is an arbitrary choice among
# degenerate states and the down-spin weight collapses.
if [ "$NCELL" -eq 3 ]; then
  NUP=$((2*CELLS*CELLS)); NDN=$((CELLS*CELLS))
  echo "    NOTE: n_cell=3 is half filling -> ferrimagnetic sector, S=$((CELLS*CELLS))/2"
else
  [ $((NELEC % 2)) -eq 0 ] || { echo "ERROR: n_cell*cells^2=$NELEC is odd, cannot set NUP=NDN"; exit 1; }
  NUP=$((NELEC/2)); NDN=$NUP
fi
ET=${ETRIAL:-$($PY -c "print(f'{-1.0*$NSITE:.6f}')")}
NFILL=$($PY -c "print(f'{$NELEC/$NSITE:.3f}')")
FOLD="L${L}n${NFILL}u${U}tA${T2}tt${T0}N${NUP}"

echo ">>> $FOLD   host=$(hostname)  cells=$CELLS  fine grid=${L}x${L}"
echo "    Lieb sites=$NSITE  k points=$((L*L))  n_cell=$NCELL  NUP=$NUP NDN=$NDN"
echo "    U(B,C)=$U  t2(B-C)=$T2  eps_A=$EPS_A  E_T=$ET  NP=$NP"
command -v mpiifort >/dev/null 2>&1 || { echo "ERROR: mpiifort not found"; exit 1; }

mkdir -p "$FOLD"
for f in *.f90 Makefile in.dat wf_lieb.py; do cp "$f" "$FOLD/"; done

sed -i.bak -E "s/parameter \(lx=[0-9]+,ly=[0-9]+/parameter (lx=${L},ly=${L}/" "$FOLD/parameter.f90"
sed -i.bak -E "s/parameter \(NUP=[0-9]+,NDN=[0-9]+/parameter (NUP=${NUP},NDN=${NDN}/" "$FOLD/parameter.f90"

$PY - "$FOLD/in.dat" <<PYEOF
import sys
fp=sys.argv[1]; L=open(fp).readlines()
L[2]=L[2].split(",")[0] + f",{$ET:.6f}\n"                              # deltau, etrial
L[4]=f"{$T0:.8f},{$T2:.8f},{$EPS_A:.8f},{0.0:.8f},-1\n"                # t0,t1,t2(=eps_A),tam
L[5]=f"{$U:.8f},{$V:.8f}\n"                                            # ud, vpd
open(fp,'w').writelines(L)
PYEOF

( cd "$FOLD"
  $PY wf_lieb.py > wf.log 2>&1 || { echo "TRIAL WF FAILED (see $FOLD/wf.log)"; tail -5 wf.log; exit 1; }
  [ "${SEED_ONLY:-0}" = "1" ] && { echo "    SEED_ONLY: trial built, not launching"; exit 0; }
  make clear >/dev/null 2>&1; make all > build.log 2>&1
  [ -f CPMC.exe ] || { echo "BUILD FAILED (see $FOLD/build.log)"; tail -6 build.log; exit 1; }
  nohup mpirun -np $NP ./CPMC.exe > nohup.out 2>&1 &
  echo "    STARTED pid=$!  -> $FOLD/nohup.out  (done when the last line is 'time_sec')"
)
