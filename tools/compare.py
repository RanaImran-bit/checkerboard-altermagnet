#!/usr/bin/env python3
"""Compare two common-schema result records (e.g. ED vs a QMC version, or two
QMC versions) observable-by-observable and report agreement.

For each scalar observable present in both records:
  delta      = value_b - value_a
  sigma      = sqrt(err_a^2 + err_b^2)              (combined stochastic error)
  z          = |delta| / sigma                       (standard deviations apart)
  verdict    = PASS if z <= ztol (default 3), else FAIL
ED records carry error 0, so z measures how many QMC sigmas the QMC value
sits from the exact answer.

Usage:
    python tools/compare.py REF.json TEST.json            # pretty table
    python tools/compare.py REF.json TEST.json --json     # machine-readable
"""
from __future__ import annotations
import argparse, json, math, sys
from pathlib import Path


def load(p):
    d = json.loads(Path(p).read_text())
    return d


def compare(ref, test, ztol=3.0):
    oa, ob = ref.get("observables", {}), test.get("observables", {})
    keys = [k for k in oa if k in ob]
    rows = []
    worst_z = 0.0
    for k in sorted(keys):
        va, ea = oa[k]["value"], (oa[k].get("error") or 0.0)
        vb, eb = ob[k]["value"], (ob[k].get("error") or 0.0)
        delta = vb - va
        sigma = math.hypot(ea, eb)
        z = abs(delta) / sigma if sigma > 0 else (0.0 if delta == 0 else math.inf)
        verdict = "PASS" if z <= ztol else "FAIL"
        worst_z = max(worst_z, z if math.isfinite(z) else 1e9)
        rows.append(dict(obs=k, ref=va, test=vb, delta=delta, sigma=sigma,
                         z=z, verdict=verdict))
    overall = "PASS" if all(r["verdict"] == "PASS" for r in rows) else "FAIL"
    return dict(ref_code=ref.get("code"), test_code=test.get("code"),
                ztol=ztol, overall=overall, worst_z=worst_z, rows=rows)


def fmt(res):
    L = []
    L.append(f"{res['ref_code']}  vs  {res['test_code']}   (ztol={res['ztol']})")
    L.append(f"{'observable':<18}{'ref':>14}{'test':>14}{'delta':>13}{'sigma':>11}{'z':>8}  verdict")
    L.append("-" * 92)
    for r in res["rows"]:
        z = "inf" if not math.isfinite(r["z"]) else f"{r['z']:.2f}"
        L.append(f"{r['obs']:<18}{r['ref']:>14.6f}{r['test']:>14.6f}"
                 f"{r['delta']:>13.6f}{r['sigma']:>11.5f}{z:>8}  {r['verdict']}")
    L.append("-" * 92)
    L.append(f"OVERALL: {res['overall']}   (worst z = {res['worst_z']:.2f})")
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("ref"); ap.add_argument("test")
    ap.add_argument("--ztol", type=float, default=3.0)
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args()
    res = compare(load(a.ref), load(a.test), a.ztol)
    print(json.dumps(res, indent=2) if a.json else fmt(res))
    return 0 if res["overall"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
