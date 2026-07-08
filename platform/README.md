# QMC Validation Platform (UI)

A brief FastAPI + React cockpit over the validation loop: compute an ED
reference, pick QMC runs, and compare energies (with error bars and
σ-distance / PASS-FAIL) — all speaking the common results schema.

## Run it

Two processes. From the repo root, in two terminals (both `source tools/env.sh`):

```bash
# 1) backend  (http://localhost:8000)
cd platform/backend
uvicorn app:app --reload --port 8000

# 2) frontend (http://localhost:5180, proxies /api -> :8000)
cd platform/frontend
npm install        # first time only
npm run dev
```

Open http://localhost:5180.

## Flow
1. **ED reference** — set lx/ly/nup/ndn/t/U → *Compute ED* (exact, via QuSpin).
2. **QMC runs** — anything under `results/<version>/<runid>/` with an `out.dat`
   (produced by `tools/run.sh`) appears in the list; click to select.
3. **Compare** — overlays ED vs the QMC run: energy bar chart with error bars
   and a per-observable table (Δ, σ, z, verdict). Total energy is the primary
   pass/fail observable (see ../docs/VALIDATION.md for why components differ).

## Features
- **Model & parameters** — single-band Hubbard or two-orbital altermagnet;
  set lattice, filling, hoppings, U/uxy/v, walkers/steps, and `bp`
  (back-propagation length; bp>0 = unbiased estimator, 0 = mixed).
- **Compute ED** + **Run QMC** live from those parameters, then **Compare**
  (energy chart with error bars + delta/sigma/z/PASS-FAIL table).
- **Parameter sweep** — sweep U/uxy/v over a list of values and overlay the
  ED curve (line) against QMC points (with error bars): E vs interaction.

## API (backend)
- `GET  /api/versions` — code versions + built status
- `GET  /api/runs`, `GET /api/runs/{code}/{run_id}` — parsed runs
- `POST /api/ed` — ED reference (hubbard or altermagnet; isolated subprocess)
- `POST /api/qmc` — run the Python CPMC port live (params incl. uxy/v/bp)
- `POST /api/sweep` — ED curve (+ optional QMC points) over one parameter
- `POST /api/compare` — diff two records
