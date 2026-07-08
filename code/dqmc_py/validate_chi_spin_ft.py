#!/usr/bin/env python3
"""Finite-T three-way ED gate for the momentum-resolved magnetic susceptibility
chi_zz(q) (work item W3, docs/PLAN_chi_spin_AHE.md):

    ED (spin_susc.ed_chi_spin, full-Fock JW Lehmann/Kubo)
    vs DQMC      (spin_susc.SpinDQMC.run_spin)
    vs CP-DQMC   (ftcpmc.FTCPMC.run_fb_stab(spin=True)), free-projection AND constrained

on the same (lx, ly, U, mu, beta, tam, t1) point, every q on the grid. Free-projection
CP-DQMC and DQMC are exact -> gated; the CONSTRAINED run reports the CP bias + <sign>
(not gated). Tolerance per q: |dev| < max(nsig * combined-stat-guess, rtol * ED).

    source tools/env.sh
    python code/dqmc_py/validate_chi_spin_ft.py --lx 2 --ly 2 --U 4 --mu 2 --beta 2 \
        --tam 0.3 --nmeas 1500 --ftnmeas 60 -o results/chi_spin_bench/ft_2x2_tam03.json
"""
from __future__ import annotations
import os, sys, argparse, json
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "ftcpqmc_py"))
from spin_susc import SpinDQMC, ed_chi_spin


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=2); ap.add_argument("--ly", type=int, default=2)
    ap.add_argument("--U", type=float, default=4.0); ap.add_argument("--mu", type=float, default=2.0)
    ap.add_argument("--beta", type=float, default=2.0); ap.add_argument("--dt", type=float, default=0.0625)
    ap.add_argument("--tam", type=float, default=0.0); ap.add_argument("--t1", type=float, default=0.0)
    ap.add_argument("--tp", type=float, default=0.0)
    ap.add_argument("--nwarm", type=int, default=300); ap.add_argument("--nmeas", type=int, default=1500)
    ap.add_argument("--ftdt", type=float, default=0.125, help="CP-DQMC Trotter step")
    ap.add_argument("--ftnw", type=int, default=20); ap.add_argument("--ftnmeas", type=int, default=40)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--rtol", type=float, default=0.10, help="rel tolerance vs ED (gate)")
    ap.add_argument("--atol", type=float, default=0.02, help="abs tolerance floor (gate)")
    ap.add_argument("--skip-cp", action="store_true", help="only ED vs DQMC")
    ap.add_argument("-o", "--out", help="write JSON record")
    a = ap.parse_args()
    if a.tp != 0.0 and not a.skip_cp:
        ap.error("--tp is not plumbed through FTCPMC; use --skip-cp with --tp")

    print(f"# chi_zz(q) finite-T gate: {a.lx}x{a.ly} U={a.U} mu={a.mu} beta={a.beta} "
          f"tam={a.tam} t1={a.t1} tp={a.tp}")

    # ---- ED (exact, full Fock JW Lehmann; zz + transverse pm) ----
    chi_ed, S_ed, dens_ed, chi_pm_ed, S_pm_ed = ed_chi_spin(
        a.lx, a.ly, a.U, a.mu, a.beta, a.tam, a.t1, a.tp, pm=True)
    print(f"  ED    density={dens_ed:.5f}")

    # ---- DQMC (exact within stats where sign ~ 1) ----
    dq = SpinDQMC(a.lx, a.ly, a.U, a.mu, a.beta, a.dt, a.tam, a.t1, a.seed, tp=a.tp)
    rd = dq.run_spin(a.nwarm, a.nmeas)
    print(f"  DQMC  density={rd['dens']:.5f}  <sign>={rd['sign']:.4f}  (NT={dq.NT})")

    rows = {"ED": chi_ed, "DQMC": np.array(rd["chi_q"])}
    S_rows = {"ED": S_ed, "DQMC": np.array(rd["S_q"])}
    pm_rows = {"ED": chi_pm_ed, "DQMC": np.array(rd["chi_pm_q"])}
    signs = {"DQMC": rd["sign"]}

    # ---- CP-DQMC: free projection (exact) + constrained (CP bias) ----
    if not a.skip_cp:
        from ftcpmc import FTCPMC
        for tag, constrained in (("CPfree", False), ("CPcons", True)):
            q = FTCPMC(a.lx, a.ly, a.U, a.mu, a.beta, a.ftdt, a.tam, a.t1, a.seed,
                       nw=a.ftnw, constrained=constrained, stab=True)
            r = q.run_fb_stab(a.ftnmeas, spin=True)
            rows[tag] = np.array(r["chi_spin_q"]); S_rows[tag] = np.array(r["S_spin_q"])
            pm_rows[tag] = np.array(r["chi_pm_q"])
            signs[tag] = r["sign"]
            print(f"  {tag} density={r['dens']:.5f}  <sign>={r['sign']:.4f}  "
                  f"npaths={r['npaths']}  (L={q.L}, dt={a.ftdt})")

    # ---- report ----
    hdr = "".join(f"{k:>10}" for k in rows)
    print(f"\n  chi_zz(q)  [per site]\n  {'q':>7} {hdr}")
    fails = []
    for kx in range(a.lx):
        for ky in range(a.ly):
            vals = "".join(f"{rows[k][kx, ky]:>10.4f}" for k in rows)
            print(f"  ({kx},{ky}) {vals}")
            ed = chi_ed[kx, ky]
            for k in rows:
                if k in ("ED", "CPcons"):        # constrained bias reported, not gated
                    continue
                dev = abs(rows[k][kx, ky] - ed)
                if dev > max(a.rtol * abs(ed), a.atol):
                    fails.append((k, kx, ky, float(dev)))
    print(f"\n  chi_pm(q)  [transverse; SU(2): = 2 chi_zz]\n  {'q':>7} {hdr}")
    for kx in range(a.lx):
        for ky in range(a.ly):
            vals = "".join(f"{pm_rows[k][kx, ky]:>10.4f}" for k in pm_rows)
            print(f"  ({kx},{ky}) {vals}")
            ed = chi_pm_ed[kx, ky]
            for k in pm_rows:
                if k in ("ED", "CPcons"):
                    continue
                dev = abs(pm_rows[k][kx, ky] - ed)
                if dev > max(a.rtol * abs(ed), a.atol):
                    fails.append((f"pm:{k}", kx, ky, float(dev)))
    print(f"\n  S^z(q)  [equal-time, per site]\n  {'q':>7} {hdr}")
    for kx in range(a.lx):
        for ky in range(a.ly):
            vals = "".join(f"{S_rows[k][kx, ky]:>10.4f}" for k in S_rows)
            print(f"  ({kx},{ky}) {vals}")
    if "CPcons" in rows:
        dev = np.abs(rows["CPcons"] - chi_ed)
        rel = dev.max() / max(np.abs(chi_ed).max(), 1e-12)
        print(f"\n  CP bias (constrained, informational): max|dev| = {dev.max():.4f} "
              f"({100*rel:.1f}% of ED max), <sign> = {signs['CPcons']:.4f}")

    verdict = "PASS" if not fails else "FAIL"
    print(f"\n{verdict}: exact methods within max({100*a.rtol:.0f}% ED, {a.atol}) at every q"
          + ("" if not fails else f" -- failures: {fails}"))
    if a.out:
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        with open(a.out, "w") as f:
            json.dump({"params": vars(a), "dens_ed": dens_ed, "signs": signs,
                       "chi": {k: v.tolist() for k, v in rows.items()},
                       "chi_pm": {k: v.tolist() for k, v in pm_rows.items()},
                       "S": {k: v.tolist() for k, v in S_rows.items()},
                       "S_pm_ed": S_pm_ed.tolist(),
                       "verdict": verdict, "fails": fails}, f, indent=1)
        print(f"wrote {a.out}")
    sys.exit(0 if verdict == "PASS" else 1)


if __name__ == "__main__":
    main()
