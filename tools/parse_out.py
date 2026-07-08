#!/usr/bin/env python3
"""Parse a CPQMC run's out.dat (+ run_meta.json) into the common results schema.

The common schema is the contract shared by every QMC version and the ED
reference, so the platform/UI can compare them directly:

    {
      "code", "run_id", "ranks", "wall_sec",     # provenance
      "model":   {nsites, nup, ndn, U, t0, ... , dt, nwalkers, bpstep},
      "observables": {                            # scalar, each value+error
         "energy_total":   {"value", "error"},
         "energy_kinetic": {"value", "error"},
         "energy_potential":{"value", "error"},
         "energy_per_site":{"value", "error"},
      },
      "kspace": [[kx,ky], ...],
      "correlations": {"Npair":[{value,error}...], "Dk":[...]},
    }

Usage:
    python tools/parse_out.py <run_dir>            # prints JSON
    python tools/parse_out.py <run_dir> -o x.json  # writes JSON
"""
from __future__ import annotations
import argparse, json, re, sys
from pathlib import Path

# out.dat label -> (schema_section, key, kind)
INT_SCALARS = {"n_site": "nsites", "n_e": "ne", "n_up": "nup", "n_dn": "ndn",
               "filling": "filling", "n_walker": "nwalkers", "bpStep": "bpstep"}
FLOAT_SCALARS = {"Ud": "U", "t0": "t0", "t1": "t1", "t2": "t2", "Vpd": "Vpd",
                 "tam": "tam", "alphat1": "alphat1", "dt": "dt"}
VAL_ERR = {"ke": "energy_kinetic", "pe": "energy_potential", "totalEn": "energy_total"}
ARRAYS = {"kSpace", "Npair", "Dk"}

_num = re.compile(r"[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?")

def _floats(line: str) -> list[float]:
    return [float(x.replace("D", "E").replace("d", "e")) for x in _num.findall(line)]

def parse_out_dat(path: Path) -> dict:
    lines = path.read_text().splitlines()
    model: dict = {}
    obs: dict = {}
    kspace: list = []
    corr: dict = {}
    i = 0
    n = len(lines)
    while i < n:
        label = lines[i].strip()
        i += 1
        if not label:
            continue
        if label in INT_SCALARS and i < n:
            v = _floats(lines[i]); i += 1
            if v: model[INT_SCALARS[label]] = int(v[0])
        elif label in FLOAT_SCALARS and i < n:
            v = _floats(lines[i]); i += 1
            if v: model[FLOAT_SCALARS[label]] = v[0]
        elif label in VAL_ERR and i < n:
            v = _floats(lines[i]); i += 1
            if len(v) >= 2:
                obs[VAL_ERR[label]] = {"value": v[0], "error": v[1]}
            elif v:
                obs[VAL_ERR[label]] = {"value": v[0], "error": None}
        elif label in ARRAYS:
            rows = []
            while i < n and lines[i].strip() and lines[i].strip() not in (
                    set(INT_SCALARS) | set(FLOAT_SCALARS) | set(VAL_ERR) | ARRAYS):
                vals = _floats(lines[i])
                if not vals:
                    break
                rows.append(vals); i += 1
            if label == "kSpace":
                kspace = [r[:2] for r in rows if len(r) >= 2]
            else:
                corr[label] = [{"value": r[0], "error": (r[1] if len(r) > 1 else None)}
                               for r in rows]
        # unknown label: skip (its data line(s) handled on next iterations)

    # derived: energy per site
    if "energy_total" in obs and model.get("nsites"):
        et = obs["energy_total"]
        ns = model["nsites"]
        obs["energy_per_site"] = {
            "value": et["value"] / ns,
            "error": (et["error"] / ns) if et["error"] is not None else None,
        }
    return {"model": model, "observables": obs, "kspace": kspace, "correlations": corr}

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("-o", "--out")
    a = ap.parse_args()
    rd = Path(a.run_dir)
    out = rd / "out.dat"
    if not out.exists():
        print(f"ERROR: {out} not found", file=sys.stderr); return 1
    rec = parse_out_dat(out)
    meta_f = rd / "run_meta.json"
    if meta_f.exists():
        rec = {**json.loads(meta_f.read_text()), **rec}
    js = json.dumps(rec, indent=2)
    if a.out:
        Path(a.out).write_text(js); print(f"wrote {a.out}")
    else:
        print(js)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
