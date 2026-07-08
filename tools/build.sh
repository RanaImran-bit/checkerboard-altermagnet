#!/usr/bin/env bash
# Build a QMC code version against the qmc-platform conda toolchain
# (gfortran + openmpi + openblas), OUT OF SOURCE.
#
#   source tools/env.sh          # once per shell (sets SDKROOT etc.)
#   tools/build.sh <version>
#
# versions:
#   benchmark-cpqmc   record/benchmark/CPQMC  (frozen reference)
#   benchmark-dqmc    record/benchmark/DQMC   (frozen reference, BSS DQMC)
#   src               code/src                (active dev tree)
#
# Frozen sources under record/ are never edited: they are copied into
# build/<version>/ and a small, documented gfortran portability shim is
# applied to the COPY (the cluster Makefiles target Intel ifort+MKL).
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VER="${1:-}"
[ -n "$VER" ] || { echo "usage: tools/build.sh <benchmark-cpqmc|benchmark-dqmc|src>"; exit 2; }

if [ -z "${CONDA_PREFIX:-}" ] || [ "$(basename "${CONDA_PREFIX:-}")" != "qmc-platform" ]; then
  echo "ERROR: activate the env first:  source tools/env.sh" >&2; exit 1
fi

FC=mpif90
# arm64 gfortran: no -mcmodel=large; relax line length; keep backtraces.
FFLAGS="-O2 -g -fbacktrace -ffree-line-length-none -fopenmp"
LDLIBS="-L${CONDA_PREFIX}/lib -lopenblas"

case "$VER" in
  benchmark-cpqmc)
    SRCDIR="$REPO/code/record/benchmark/CPQMC"
    ORDER="parameter jiekou mc2duph libuph cpPara cpOut cp"
    EXE="CPMC.exe" ;;
  benchmark-dqmc)
    SRCDIR="$REPO/code/record/benchmark/DQMC"
    ORDER="BSS BSSOutput AppBSS"
    EXE="SOCh.exe" ;;
  src|src-small)
    SRCDIR="$REPO/code/src"
    ORDER="cpParameter z_libs z_jiekou t Vee Vph z_overlap z_stabilize z_walkers z_energy meas cpOut cpParallel cpCore cpMain"
    EXE="CPMC.exe" ;;
  *) echo "unknown version: $VER" >&2; exit 2 ;;
esac

BD="$REPO/build/$VER"
rm -rf "$BD"; mkdir -p "$BD"
cp "$SRCDIR"/*.f90 "$BD"/

# ---- gfortran portability shims (applied to staged copies only) ----
apply_shims() {
  case "$VER" in
    benchmark-cpqmc)
      # ifort auto-pads character-array constructors; gfortran needs a typed
      # constructor. kNames(...) = (/ 'gx_up','n_up',... /) -> typed form.
      sed -i '' 's#= (/ \&#= (/ character(len=32) :: \&#' "$BD/mc2duph.f90" ;;
    src-small)
      # small cluster for ED comparison. Override via env (defaults: 2x2, 4+4,
      # 1000 walkers):
      #   QMC_LX QMC_LY QMC_NUP QMC_NDN QMC_NWLKRS
      # (QMC_NWLKRS raises walker count for tighter error bars on ED checks.)
      LXv="${QMC_LX:-2}"; LYv="${QMC_LY:-2}"; NUPv="${QMC_NUP:-4}"; NDNv="${QMC_NDN:-4}"
      NWv="${QMC_NWLKRS:-1000}"
      sed -i '' "s#lx=4, ly=4#lx=${LXv}, ly=${LYv}#" "$BD/cpParameter.f90"
      sed -i '' "s#NUP=16, NDN=16#NUP=${NUPv}, NDN=${NDNv}#" "$BD/cpParameter.f90"
      sed -i '' "s#NWLKRS=1000#NWLKRS=${NWv}#" "$BD/cpParameter.f90"
      echo "[build:src-small] patched cpParameter: lx=${LXv} ly=${LYv}, NUP=${NUPv} NDN=${NDNv}, NWLKRS=${NWv}" ;;
  esac
}
apply_shims

cd "$BD"
echo "[build:$VER] FC=$FC  flags=$FFLAGS"
for f in $ORDER; do
  printf '  cc %-14s ' "$f.f90"
  if $FC $FFLAGS -c "$f.f90" -o "$f.o" > "$f.log" 2>&1; then echo ok
  else echo FAIL; echo "---- $f errors ----"; grep -iE 'error' "$f.log" | head -20; exit 1; fi
done

OBJS=""; for f in $ORDER; do OBJS="$OBJS $f.o"; done
echo "[build:$VER] link -> $EXE"
$FC $FFLAGS -o "$EXE" $OBJS $LDLIBS 2> link.log || { grep -v 'duplicate -rpath' link.log; exit 1; }
echo "[build:$VER] OK -> $BD/$EXE"
