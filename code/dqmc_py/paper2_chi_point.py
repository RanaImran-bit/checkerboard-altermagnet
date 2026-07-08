#!/usr/bin/env python3
"""Paper-2 demo driver: momentum-resolved magnetic susceptibilities chi_zz(q) and
chi_pm(q) for the combined altermagnet model (spin-dependent NN tA + NNN t')
at ONE (n, tA, t') point with ONE method. The model == am_hopping(t=1, tam=tA, t1=t'):

    eps_up(k) = -2[(t-tA)cos kx + (t+tA)cos ky] + 4 t' sin kx sin ky   (paper Eq. 9)
    eps_dn(k) = eps_up(R90 k)                                          (paper Eq. 10)

Methods:
  dqmc    finite-T BSS (exact; sign-limited away from half filling -- sign reported)
  cpdqmc  finite-T constrained-path (sign-controlled workhorse)
  cpqmc   T=0 CPQMC, nearest CLOSED-SHELL filling (windowed one-sided Kubo;
          x2 to compare with the finite-T convention)

mu for the finite-T methods: Hartree estimate mu0(n) + U n/2, refined by one
secant step on short runs. Writes one JSON per (point, method).

    python code/dqmc_py/paper2_chi_point.py --method cpdqmc --n 0.9 --ta 0.5 --tt 0.1 \
        --lx 8 --ly 8 --beta 5 -o results/paper2_demo/cpdqmc_n0.9_ta0.5_tt0.1.json
"""
from __future__ import annotations
import os, sys, argparse, json
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "ftcpqmc_py"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "pyqmc"))
from cpqmc import am_hopping


def free_levels(lx, ly, ta, tt):
    """Sorted single-particle levels of one spin species (both species share the
    spectrum: eps_dn(k) = eps_up(R90 k))."""
    Ku, _ = am_hopping(lx, ly, 1.0, ta, tt)
    return np.sort(np.linalg.eigvalsh(Ku))


def mu_hartree(lx, ly, ta, tt, n, U, beta):
    """mu estimate: noninteracting mu at filling n (fine k-grid, T=1/beta) + U n/2."""
    M = 512
    k = 2 * np.pi * np.arange(M) / M
    kx, ky = np.meshgrid(k, k, indexing="ij")
    eu = -2 * ((1 - ta) * np.cos(kx) + (1 + ta) * np.cos(ky)) + 4 * tt * np.sin(kx) * np.sin(ky)
    ed = -2 * ((1 + ta) * np.cos(kx) + (1 - ta) * np.cos(ky)) - 4 * tt * np.sin(kx) * np.sin(ky)
    T = 1.0 / beta
    lo, hi = eu.min() - 5, eu.max() + 5
    for _ in range(60):
        mu = 0.5 * (lo + hi)
        f = lambda e: 1.0 / (np.exp(np.clip((e - mu) / T, -40, 40)) + 1.0)
        nn = f(eu).mean() + f(ed).mean()
        lo, hi = (mu, hi) if nn < n else (lo, mu)
    return 0.5 * (lo + hi) + 0.5 * U * n


def closed_shell_filling(lx, ly, ta, tt, n_target, min_gap=1e-6):
    """Nearest filling to n_target with a nondegenerate Fermi level (closed shell)."""
    w = free_levels(lx, ly, ta, tt)
    nsite = lx * ly
    n0 = int(round(0.5 * n_target * nsite))
    best = None
    for dn in range(0, nsite // 4):
        for cand in ({n0 + dn, n0 - dn}):
            if cand < 1 or cand >= nsite:
                continue
            gap = w[cand] - w[cand - 1]
            if gap > min_gap:
                sc = (abs(cand - n0), -gap)
                if best is None or sc < best[0]:
                    best = (sc, cand, gap)
        if best is not None and abs(best[1] - n0) <= dn:
            break
    _, nup, gap = best
    return nup, 2.0 * nup / nsite, float(gap)


def run_dqmc(a, mu):
    from spin_susc import SpinDQMC
    q = SpinDQMC(a.lx, a.ly, a.U, mu, a.beta, a.dt, a.ta, a.tt, a.seed)
    r = q.run_spin(a.nwarm, a.nmeas)
    return dict(chi_q=np.real(r["chi_q"]).tolist(), chi_pm_q=np.real(r["chi_pm_q"]).tolist(),
                S_q=np.real(r["S_q"]).tolist(), S_pm_q=np.real(r["S_pm_q"]).tolist(),
                dens=r["dens"], sign=r["sign"], beta=a.beta, mu=mu)


def run_cpdqmc(a, mu, nmeas=None):
    from ftcpmc import FTCPMC
    q = FTCPMC(a.lx, a.ly, a.U, mu, a.beta, a.ftdt, a.ta, a.tt, a.seed,
               nw=a.ftnw, constrained=True, stab=True)
    r = q.run_fb_stab(nmeas or a.ftnmeas, spin=(nmeas is None))
    out = dict(dens=r["dens"], sign=r["sign"], npaths=r["npaths"], beta=a.beta, mu=mu)
    if nmeas is None:
        out.update(chi_q=r["chi_spin_q"], chi_pm_q=r["chi_pm_q"],
                   S_q=r["S_spin_q"], S_pm_q=r["S_pm_q"])
    return out


def run_cpqmc(a):
    from cpqmc import CPMC
    nup, n_act, gap = closed_shell_filling(a.lx, a.ly, a.ta, a.tt, a.n)
    Ku, Kd = am_hopping(a.lx, a.ly, 1.0, a.ta, a.tt)
    q = CPMC(a.lx, a.ly, nup, nup, t=1.0, U=a.U, dt=a.t0dt, nwalkers=a.nw,
             seed=a.seed, K=Ku, K_dn=Kd)
    r = q.run_bp_chi_spin(nequil=a.nequil, nblocks=a.nblocks, bp=a.bp)
    # x2: one-sided T=0 windowed Kubo -> finite-T convention (docs/PLAN_chi_spin_AHE.md)
    return dict(chi_q=(2 * np.array(r["chi_q"])).tolist(),
                chi_pm_q=(2 * np.array(r["chi_pm_q"])).tolist(),
                chi_q_err=(2 * np.array(r["chi_q_err"])).tolist(),
                S_q=np.array(r["Ctau_q"])[0].tolist(),        # C_q(tau=0) = S^z(q)
                S_pm_q=np.array(r["Ctau_pm_q"])[0].tolist(),
                dens=n_act, nup=nup, shell_gap=gap, bp=a.bp, dt=a.t0dt,
                note="T=0 windowed chi x2 (finite-T convention); closed-shell filling")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--method", choices=["dqmc", "cpdqmc", "cpqmc"], required=True)
    ap.add_argument("--n", type=float, required=True); ap.add_argument("--ta", type=float, required=True)
    ap.add_argument("--tt", type=float, required=True)
    ap.add_argument("--lx", type=int, default=8); ap.add_argument("--ly", type=int, default=8)
    ap.add_argument("--U", type=float, default=4.0); ap.add_argument("--beta", type=float, default=5.0)
    ap.add_argument("--dt", type=float, default=0.125)      # DQMC Trotter
    ap.add_argument("--ftdt", type=float, default=0.125)    # CP-DQMC Trotter
    ap.add_argument("--t0dt", type=float, default=0.02)     # CPQMC Trotter
    ap.add_argument("--nwarm", type=int, default=200); ap.add_argument("--nmeas", type=int, default=600)
    ap.add_argument("--ftnw", type=int, default=20); ap.add_argument("--ftnmeas", type=int, default=60)
    ap.add_argument("--nw", type=int, default=300); ap.add_argument("--nequil", type=int, default=300)
    ap.add_argument("--nblocks", type=int, default=20); ap.add_argument("--bp", type=int, default=60)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()

    rec = dict(params=vars(a))
    if a.method == "cpqmc":
        rec.update(run_cpqmc(a))
    else:
        mu = mu_hartree(a.lx, a.ly, a.ta, a.tt, a.n, a.U, a.beta)
        # one secant refinement on a short sign-controlled run
        r1 = run_cpdqmc(a, mu, nmeas=6)
        dmu = 0.4
        r2 = run_cpdqmc(a, mu + dmu, nmeas=6)
        if abs(r2["dens"] - r1["dens"]) > 1e-4:
            mu = mu + dmu * (a.n - r1["dens"]) / (r2["dens"] - r1["dens"])
        rec["mu_tuning"] = dict(mu_hartree=r1["mu"], n1=r1["dens"], n2=r2["dens"], mu_final=mu)
        rec.update(run_dqmc(a, mu) if a.method == "dqmc" else run_cpdqmc(a, mu))
    # peak diagnostics (both channels)
    fold = lambda v: v - 2.0 if v > 1.0 else v          # q/pi in (-1, 1]
    for ch in ("chi_q", "chi_pm_q"):
        g = np.array(rec[ch])
        kx, ky = np.unravel_index(np.argmax(g), g.shape)
        rec[ch + "_peak"] = dict(qx_over_pi=fold(2.0 * kx / a.lx),
                                 qy_over_pi=fold(2.0 * ky / a.ly),
                                 value=float(g[kx, ky]))
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w") as f:
        json.dump(rec, f, indent=1)
    pz = rec["chi_q_peak"]; pp = rec["chi_pm_q_peak"]
    print(f"{a.method} n={a.n} tA={a.ta} t'={a.tt}: dens={rec.get('dens'):.4f} "
          f"sign={rec.get('sign', 1):.3f}  chi_zz peak {pz['value']:.4f} at "
          f"({pz['qx_over_pi']:.2f},{pz['qy_over_pi']:.2f})pi   chi_pm peak {pp['value']:.4f} at "
          f"({pp['qx_over_pi']:.2f},{pp['qy_over_pi']:.2f})pi -> wrote {a.out}")


if __name__ == "__main__":
    main()
