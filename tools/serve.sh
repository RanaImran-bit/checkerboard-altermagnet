#!/usr/bin/env bash
# Launch the validation-platform backend (FastAPI :8000) and frontend
# (Vite :5180) together. Ctrl-C stops both.
#
#   source tools/env.sh
#   tools/serve.sh
#   # then open http://localhost:5180
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
[ "$(basename "${CONDA_PREFIX:-}")" = "qmc-platform" ] || { echo "source tools/env.sh first"; exit 1; }

# single-thread keeps the ED subprocess clear of the openblas/quspin OpenMP clash
export OMP_NUM_THREADS=1 KMP_DUPLICATE_LIB_OK=TRUE

cleanup() { kill "${BACK_PID:-}" "${FRONT_PID:-}" 2>/dev/null || true; }
trap cleanup EXIT INT TERM

( cd "$REPO/platform/backend" && uvicorn app:app --host 127.0.0.1 --port 8000 ) &
BACK_PID=$!
[ -d "$REPO/platform/frontend/node_modules" ] || ( cd "$REPO/platform/frontend" && npm install )
( cd "$REPO/platform/frontend" && npm run dev -- --host ) &
FRONT_PID=$!

echo "backend  http://localhost:8000   (pid $BACK_PID)"
echo "frontend http://localhost:5180   (pid $FRONT_PID)"
wait
