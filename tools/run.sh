#!/usr/bin/env bash
# Run a built QMC version and collect outputs into results/<version>/<run-id>/.
#
#   source tools/env.sh
#   tools/build.sh <version>
#   tools/run.sh <version> [in.dat] [nranks]
#
# Defaults: the version's shipped in.dat, 2 MPI ranks.
# Sets the macOS stack limit high enough for the code's large arrays
# (the cluster build relied on -mcmodel=large; locally we lift ulimit -s).
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VER="${1:-}"
[ -n "$VER" ] || { echo "usage: tools/run.sh <version> [in.dat] [nranks]"; exit 2; }
INPUT="${2:-}"
NRANKS="${3:-2}"

[ "$(basename "${CONDA_PREFIX:-}")" = "qmc-platform" ] || { echo "ERROR: source tools/env.sh first" >&2; exit 1; }

BD="$REPO/build/$VER"
EXE="$(ls "$BD"/*.exe 2>/dev/null | head -1)"
[ -x "$EXE" ] || { echo "ERROR: no executable in $BD — run tools/build.sh $VER" >&2; exit 1; }

# default input = the version's source in.dat
case "$VER" in
  benchmark-cpqmc) SRCDIR="$REPO/code/record/benchmark/CPQMC" ;;
  benchmark-dqmc)  SRCDIR="$REPO/code/record/benchmark/DQMC" ;;
  src)             SRCDIR="$REPO/code/src" ;;
esac
[ -n "$INPUT" ] || INPUT="$SRCDIR/in.dat"
[ -f "$INPUT" ] || { echo "ERROR: input not found: $INPUT" >&2; exit 1; }

RUNID="$(date +%Y%m%d-%H%M%S)"
RD="$REPO/results/$VER/$RUNID"
mkdir -p "$RD"
cp "$EXE" "$RD/"; cp "$INPUT" "$RD/in.dat"

# large-array / OpenMP stack headroom
ulimit -s 65520 2>/dev/null || true
export OMP_STACKSIZE=512M GOMP_STACKSIZE=512000

cd "$RD"
echo "[run:$VER] id=$RUNID ranks=$NRANKS input=$(basename "$INPUT")"
START=$(python -c 'import time;print(time.time())')
mpirun --oversubscribe -np "$NRANKS" "./$(basename "$EXE")" > stdout.log 2>&1 || true
END=$(python -c 'import time;print(time.time())')
WALL=$(python -c "print(f'{$END-$START:.2f}')")

# record run metadata
cat > run_meta.json <<EOF
{"code":"$VER","run_id":"$RUNID","ranks":$NRANKS,"wall_sec":$WALL,
 "input":"$(basename "$INPUT")","exe":"$(basename "$EXE")"}
EOF

echo "[run:$VER] done in ${WALL}s -> $RD"
ls -1 "$RD"
[ -f out.dat ] && echo "[run:$VER] out.dat produced ($(wc -l < out.dat) lines)" \
              || echo "[run:$VER] WARNING: no out.dat (check stdout.log)"
