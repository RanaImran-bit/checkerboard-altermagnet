#!/usr/bin/env python3
"""FastAPI backend for the QMC validation platform.

Serves the common-schema results loop to the UI:
  - list code versions and whether they are built
  - list/parse existing runs under results/
  - compute an ED reference (subprocess -> ed/hubbard_ed.py, isolated so the
    server process never imports quspin / hits the OpenMP clash)
  - compare any two common-schema records (delta / sigma / z / PASS-FAIL)

Run:
  source tools/env.sh
  uvicorn app:app --reload --port 8000      # from platform/backend/
"""
from __future__ import annotations
import json
import subprocess
import sys
from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

REPO = Path(__file__).resolve().parents[2]
TOOLS = REPO / "tools"
RESULTS = REPO / "results"
BUILD = REPO / "build"
sys.path.insert(0, str(TOOLS))
import parse_out          # noqa: E402  (common-schema parser)
import compare as cmp     # noqa: E402  (diff engine)

VERSIONS = {
    "benchmark-cpqmc": REPO / "code/record/benchmark/CPQMC",
    "benchmark-dqmc": REPO / "code/record/benchmark/DQMC",
    "src": REPO / "code/src",
}

app = FastAPI(title="QMC Validation Platform")
app.add_middleware(
    CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"]
)


@app.get("/api/health")
def health():
    return {"ok": True, "repo": str(REPO)}


@app.get("/api/versions")
def versions():
    out = []
    for vid, src in VERSIONS.items():
        exe = list((BUILD / vid).glob("*.exe")) if (BUILD / vid).exists() else []
        out.append({"id": vid, "src": str(src.relative_to(REPO)), "built": bool(exe)})
    return out


def _parse_run(code: str, run_dir: Path) -> dict | None:
    out = run_dir / "out.dat"
    if not out.exists():
        return None
    rec = parse_out.parse_out_dat(out)
    meta = run_dir / "run_meta.json"
    if meta.exists():
        rec = {**json.loads(meta.read_text()), **rec}
    rec.setdefault("code", code)
    rec.setdefault("run_id", run_dir.name)
    return rec


@app.get("/api/runs")
def runs():
    """All parsed runs across versions (summary)."""
    out = []
    if not RESULTS.exists():
        return out
    for code_dir in sorted(RESULTS.iterdir()):
        if not code_dir.is_dir():
            continue
        for run_dir in sorted(code_dir.iterdir(), reverse=True):
            rec = _parse_run(code_dir.name, run_dir)
            if rec:
                et = rec.get("observables", {}).get("energy_total", {})
                out.append({
                    "code": rec["code"], "run_id": rec["run_id"],
                    "model": rec.get("model", {}),
                    "energy_total": et.get("value"), "error": et.get("error"),
                    "wall_sec": rec.get("wall_sec"),
                })
    return out


@app.get("/api/runs/{code}/{run_id}")
def run_detail(code: str, run_id: str):
    rec = _parse_run(code, RESULTS / code / run_id)
    if not rec:
        raise HTTPException(404, "run or out.dat not found")
    return rec


class EDReq(BaseModel):
    model: str = "hubbard"          # "hubbard" | "altermagnet"
    lx: int = 4
    ly: int = 4
    nup: int = 1
    ndn: int = 1
    t: float = 1.0
    U: float = 3.0
    t1: float = -1.0
    t2: float = -1.0
    t3: float = -1.0
    t4: float = -1.0
    uxy: float = 0.0
    v: float = 0.0


def _ed_record(p: dict) -> dict:
    """Run ED for a param dict; dispatch hubbard vs altermagnet."""
    if p.get("model") == "altermagnet":
        cmd = [sys.executable, str(REPO / "ed/altermagnet_ed.py"),
               "--lx", str(p["lx"]), "--ly", str(p["ly"]),
               "--nup", str(p["nup"]), "--ndn", str(p["ndn"]),
               "--t1", str(p["t1"]), "--t2", str(p["t2"]),
               "--t3", str(p["t3"]), "--t4", str(p["t4"]),
               "--uxx", str(p["U"]), "--uxy", str(p["uxy"]), "--v", str(p["v"])]
    else:
        cmd = [sys.executable, str(REPO / "ed/hubbard_ed.py"),
               "--lx", str(p["lx"]), "--ly", str(p["ly"]),
               "--nup", str(p["nup"]), "--ndn", str(p["ndn"]),
               "--t", str(p.get("t", 1.0)), "--U", str(p["U"])]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    if r.returncode != 0:
        raise HTTPException(500, f"ED failed: {r.stderr[-500:]}")
    return json.loads(r.stdout)


@app.post("/api/ed")
def ed(req: EDReq):
    """Exact-diagonalization reference (isolated subprocess)."""
    return _ed_record(req.model_dump())


class CompareReq(BaseModel):
    ref: dict
    test: dict
    ztol: float = 3.0


@app.post("/api/compare")
def compare(req: CompareReq):
    return cmp.compare(req.ref, req.test, req.ztol)


class QMCReq(BaseModel):
    model: str = "hubbard"          # "hubbard" | "altermagnet"
    lx: int = 4
    ly: int = 4
    nup: int = 1
    ndn: int = 1
    U: float = 3.0
    t: float = 1.0
    t1: float = -1.0
    t2: float = -1.0
    t3: float = -1.0
    t4: float = -1.0
    uxy: float = 0.0
    v: float = 0.0
    dt: float = 0.01
    nwalkers: int = 200
    nequil: int = 150
    nmeas: int = 300
    bp: int = 0                     # back-propagation length (0 = mixed estimator)


def _qmc_record(p: dict, tag: str = "_tmp") -> dict:
    """Run the Python CPMC port for a param dict; returns the common schema."""
    out = REPO / "results" / f"_qmc_{tag}.json"
    cmd = [sys.executable, str(REPO / "pyqmc/cpqmc.py"),
           "--model", p["model"], "--lx", str(p["lx"]), "--ly", str(p["ly"]),
           "--nup", str(p["nup"]), "--ndn", str(p["ndn"]), "--U", str(p["U"]),
           "--t", str(p.get("t", 1.0)), "--t1", str(p["t1"]), "--t2", str(p["t2"]),
           "--t3", str(p["t3"]), "--t4", str(p["t4"]),
           "--uxy", str(p["uxy"]), "--v", str(p["v"]), "--dt", str(p["dt"]),
           "--nw", str(p["nwalkers"]), "--nequil", str(p["nequil"]),
           "--nmeas", str(p["nmeas"]), "--bp", str(p.get("bp", 0)), "-o", str(out)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    if r.returncode != 0 or not out.exists():
        raise HTTPException(500, f"QMC failed: {r.stderr[-600:] or r.stdout[-600:]}")
    rec = json.loads(out.read_text())
    out.unlink(missing_ok=True)
    rec["estimator"] = f"back-prop (bp={p['bp']})" if p.get("bp", 0) > 0 else "mixed"
    return rec


@app.post("/api/qmc")
def qmc(req: QMCReq):
    """Run the (validated) Python CPMC port on demand. Set bp>0 for the
    unbiased back-propagated estimator."""
    rec = _qmc_record(req.model_dump())
    if req.model == "altermagnet" and req.v != 0 and req.bp == 0:
        rec["note"] = ("Mixed estimator: the neighbour-v energy carries a bias "
                       "that grows with v. Enable back-propagation (bp>0) for an "
                       "unbiased result.")
    return rec


class CorrReq(BaseModel):
    lx: int = 2
    ly: int = 2
    nup: int = 4
    ndn: int = 4
    U: float = 2.0
    t1: float = -1.0
    t2: float = -1.0
    t3: float = -1.0
    t4: float = -1.0
    uxy: float = 0.0
    v: float = 0.0
    bp: int = 12
    nmeas: int = 360


@app.post("/api/correlations")
def correlations(req: CorrReq):
    """Compare the equal-time Green's function and spin/charge correlations
    between ED and the (back-propagated) Python CPQMC, element by element."""
    import numpy as np
    p = req.model_dump()
    ed_out = REPO / "results" / "_corr_ed.json"
    qm_out = REPO / "results" / "_corr_qmc.json"
    ed_cmd = [sys.executable, str(REPO / "ed/altermagnet_ed.py"), "--corr",
              "--lx", str(p["lx"]), "--ly", str(p["ly"]), "--nup", str(p["nup"]),
              "--ndn", str(p["ndn"]), "--uxx", str(p["U"]), "--uxy", str(p["uxy"]),
              "--v", str(p["v"]), "--t1", str(p["t1"]), "--t2", str(p["t2"]),
              "--t3", str(p["t3"]), "--t4", str(p["t4"]), "-o", str(ed_out)]
    qm_cmd = [sys.executable, str(REPO / "pyqmc/cpqmc.py"), "--corr", "--model", "altermagnet",
              "--lx", str(p["lx"]), "--ly", str(p["ly"]), "--nup", str(p["nup"]),
              "--ndn", str(p["ndn"]), "--U", str(p["U"]), "--uxy", str(p["uxy"]),
              "--v", str(p["v"]), "--t1", str(p["t1"]), "--t2", str(p["t2"]),
              "--t3", str(p["t3"]), "--t4", str(p["t4"]), "--bp", str(p["bp"]),
              "--nmeas", str(p["nmeas"]), "-o", str(qm_out)]
    for cmd, name in ((ed_cmd, "ED"), (qm_cmd, "QMC")):
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
        if r.returncode != 0:
            raise HTTPException(500, f"{name} correlations failed: {r.stderr[-500:]}")
    edr = json.loads(ed_out.read_text()); qmr = json.loads(qm_out.read_text())
    summary = []
    for key, label in (("green_up", "Green G↑"), ("nn", "charge ⟨nᵢnⱼ⟩"),
                       ("szsz", "spin ⟨SᶻᵢSᶻⱼ⟩")):
        E = np.array(edr[key]); Q = np.array(qmr[key]); d = np.abs(E - Q)
        summary.append({"obs": label, "max_dev": float(d.max()), "mean_dev": float(d.mean()),
                        "ed_range": [float(E.min()), float(E.max())],
                        "verdict": "PASS" if d.max() < 0.05 else "FAIL"})
    for f in (ed_out, qm_out):
        f.unlink(missing_ok=True)
    return {"nsites": edr["nsites"], "summary": summary,
            "ed": {k: edr[k] for k in ("green_up", "nn", "szsz")},
            "qmc": {k: qmr[k] for k in ("green_up", "nn", "szsz")}}


class BenchReq(BaseModel):
    lx: int = 2
    ly: int = 2
    nup: int = 4
    ndn: int = 4
    uxx: float = 2.0
    fort_nblk: int = 2
    fort_nblkstps: int = 20
    py_nw: int = 100
    py_nsteps: int = 20
    ranks: int = 2


@app.post("/api/benchmark")
def benchmark(req: BenchReq):
    """Fortran vs Python CPQMC speed (wall time + throughput)."""
    p = req.model_dump()
    out = REPO / "results" / "_bench.json"
    cmd = [sys.executable, str(REPO / "tools/benchmark.py"),
           "--lx", str(p["lx"]), "--ly", str(p["ly"]), "--nup", str(p["nup"]),
           "--ndn", str(p["ndn"]), "--uxx", str(p["uxx"]),
           "--fort-nblk", str(p["fort_nblk"]), "--fort-nblkstps", str(p["fort_nblkstps"]),
           "--py-nw", str(p["py_nw"]), "--py-nsteps", str(p["py_nsteps"]),
           "--ranks", str(p["ranks"]), "-o", str(out)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    if r.returncode != 0 or not out.exists():
        raise HTTPException(500, f"benchmark failed: {r.stderr[-600:] or r.stdout[-600:]}")
    res = json.loads(out.read_text()); out.unlink(missing_ok=True)
    return res


class SweepReq(BaseModel):
    base: QMCReq                    # base parameters
    var: str = "U"                  # swept variable: "U" | "uxy" | "v"
    values: list[float] = [0.0, 0.5, 1.0, 1.5, 2.0]
    run_qmc: bool = True            # also run QMC at each point (slower)


@app.post("/api/sweep")
def sweep(req: SweepReq):
    """Sweep one interaction parameter and return the ED curve (always) plus
    optional QMC points, so the UI can overlay E vs parameter."""
    pts = []
    for i, x in enumerate(req.values):
        p = req.base.model_dump()
        p[req.var] = x
        ed_rec = _ed_record(p)
        row = {"x": x, "ed": ed_rec["observables"]["energy_total"]["value"]}
        if req.run_qmc:
            q = _qmc_record(p, tag=f"sweep{i}")
            et = q["observables"]["energy_total"]
            row["qmc"] = et["value"]; row["qmc_err"] = et.get("error")
        pts.append(row)
    return {"var": req.var, "points": pts}


class PairingReq(BaseModel):
    """d-wave pairing-vertex benchmark vs full-Fock ED. Defaults to the
    NON-DEGENERATE 4x2 / 3+3 case (6 electrons on 8 sites), where the
    free-electron-trial CP-AFQMC reproduces BOTH the energy and the equal-time
    d-wave vertex of ED -- unlike the frustrated half-filled (degenerate) point."""
    lx: int = 4
    ly: int = 2
    nup: int = 3
    ndn: int = 3
    U: float = 4.0
    t0: float = 1.0
    tam: float = 0.2
    t1: float = 0.3
    etas: str = "0.0,0.5"          # 0 = free bra; >0 = back-propagated AGP/BCS bra
    dt: float = 0.02
    bp: int = 28
    nblocks: int = 30
    npop: int = 12
    nw: int = 60
    seed: int = 1


@app.post("/api/pairing")
def pairing(req: PairingReq):
    """Benchmark the back-propagated AGP/BCS d-wave pairing VERTEX (and the energy)
    against full-Fock ED on a non-degenerate cluster. Runs pyqmc/agp_bp_vertex.py
    in an isolated subprocess (it computes the ED reference internally via quspin),
    then scores energy (z) and the d-wave vertex (fraction of ED) per eta."""
    p = req.model_dump()
    out = REPO / "results" / "_pairing.json"
    cmd = [sys.executable, str(REPO / "pyqmc/agp_bp_vertex.py"),
           "--lx", str(p["lx"]), "--ly", str(p["ly"]), "--nup", str(p["nup"]),
           "--ndn", str(p["ndn"]), "--U", str(p["U"]), "--t0", str(p["t0"]),
           "--tam", str(p["tam"]), "--t1", str(p["t1"]), "--etas", p["etas"],
           "--dt", str(p["dt"]), "--bp", str(p["bp"]), "--nblocks", str(p["nblocks"]),
           "--npop", str(p["npop"]), "--nw", str(p["nw"]), "--seed", str(p["seed"]),
           "--ed", "-o", str(out)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=1200)
    if r.returncode != 0 or not out.exists():
        raise HTTPException(500, f"pairing benchmark failed: {r.stderr[-700:] or r.stdout[-700:]}")
    rec = json.loads(out.read_text())
    out.unlink(missing_ok=True)
    ed = rec.get("ed") or {}
    e_ed = ed.get("energy"); v_ed = ed.get("S_vertex")
    rows = []
    for row in rec["rows"]:
        z = (abs(row["energy"] - e_ed) / row["energy_err"]
             if (e_ed is not None and row["energy_err"] > 0) else None)
        frac = (row["vertex"] / v_ed if (v_ed not in (None, 0)) else None)
        rows.append({
            "eta": row["eta"], "bra": "free" if row["eta"] == 0 else f"AGP η={row['eta']}",
            "energy": row["energy"], "energy_err": row["energy_err"],
            "vertex": row["vertex"], "vertex_err": row["vertex_err"],
            "z": z, "energy_pass": (z is not None and z < 3.0),
            "vertex_frac": frac,
        })
    return {"params": {k: p[k] for k in ("lx", "ly", "nup", "ndn", "U", "tam", "t1", "dt", "bp")},
            "degenerate": (p["nup"] == p["ndn"] and 2 * p["nup"] == p["lx"] * p["ly"]),
            "ed": {"energy": e_ed, "S_full": ed.get("S_full"), "S_vertex": v_ed},
            "sign": rec.get("sign"), "rows": rows}
