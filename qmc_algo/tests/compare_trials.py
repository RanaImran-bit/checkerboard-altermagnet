#!/usr/bin/env python3
"""Head-to-head: free determinant vs ADAPTIVE vs CASCI vs EVOLUTIONARY trial.

Same lattice, same U, same orbital basis W (K-eigenbasis = occupation order), same
CPMC settings, same walker re-initialization. We report the constrained-path bias
E_CPMC - E_ED for each trial-construction method as a function of the number of
determinants it uses, so the comparison is at MATCHED cost.

The question this answers: does evolutionary selected-CI beat adaptive (a single
self-consistent determinant) and CASCI (all configs in a natural-orbital active
window) -- the two existing ED-free trial builders in this repo?

    python compare_trials.py --lx 2 --ly 2 --nup 2 --ndn 2 --U 4
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np

_HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(_HERE, "..", "src"))
sys.path.insert(0, os.path.join(_HERE, "..", "..", "pyqmc"))
sys.path.insert(0, os.path.join(_HERE, "..", "..", "ed"))

from evo_trial import evolve_trial, evolved_multidet, selected_ci_solve  # noqa
from casci import transform_hubbard, casci_trial
from cpqmc import CPMC, square_hopping


def run_multidet(lx, ly, nup, ndn, t, U, dt, nw, nequil, nmeas, seed,
                 ups, dns, coef):
    """CPMC with an externally-built multideterminant trial, walkers re-initialized
    on the leading determinant (finite starting overlap)."""
    q = CPMC(lx, ly, nup, ndn, t=t, U=U, dt=dt, nwalkers=nw, seed=seed,
             trial="multidet")
    q.trial.set_multidet(ups, dns, coef)
    lu = np.asarray(ups[0], float); ld = np.asarray(dns[0], float)
    q.walkers.phi_up = np.stack([lu.copy() for _ in range(q.walkers.nw)])
    q.walkers.phi_dn = np.stack([ld.copy() for _ in range(q.walkers.nw)])
    q.walkers.w = np.ones(q.walkers.nw)
    q.trial._refresh_olp(q.walkers)
    if q.walkers.olp[0] < 0:
        q.trial.coef = -q.trial.coef; q.trial._refresh_olp(q.walkers)
    return q.run(nequil=nequil, nmeas=nmeas)


def run_single(lx, ly, nup, ndn, t, U, dt, nw, nequil, nmeas, seed, mode,
               trial_mix=0.5, trial_every=10):
    q = CPMC(lx, ly, nup, ndn, t=t, U=U, dt=dt, nwalkers=nw, seed=seed,
             trial=mode, trial_mix=trial_mix, trial_every=trial_every)
    return q.run(nequil=nequil, nmeas=nmeas)


def ed_energy(lx, ly, nup, ndn, t, U):
    K = square_hopping(lx, ly, t); n = lx * ly
    _, W = np.linalg.eigh(K)
    from itertools import combinations
    full = [(u, d) for u in combinations(range(n), nup)
                    for d in combinations(range(n), ndn)]
    h1, Up = transform_hubbard(K, U, W)
    E0, _, _ = selected_ci_solve(full, h1, Up)
    try:
        from hubbard_ed import build
        _, H, _, _ = build(lx, ly, nup, ndn, t, U)
        E0 = float(np.linalg.eigvalsh(H.toarray())[0])
    except Exception:
        pass
    return K, W, h1, Up, E0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=2); ap.add_argument("--ly", type=int, default=2)
    ap.add_argument("--nup", type=int, default=2); ap.add_argument("--ndn", type=int, default=2)
    ap.add_argument("--U", type=float, default=4.0); ap.add_argument("--t", type=float, default=1.0)
    ap.add_argument("--dt", type=float, default=0.01); ap.add_argument("--nw", type=int, default=400)
    ap.add_argument("--nequil", type=int, default=150); ap.add_argument("--nmeas", type=int, default=300)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    n = a.lx * a.ly
    K, W, h1, Up, E0 = ed_energy(a.lx, a.ly, a.nup, a.ndn, a.t, a.U)
    print(f"# {a.lx}x{a.ly}, nup={a.nup} ndn={a.ndn}, U={a.U}, t={a.t}; ED E0 = {E0:.6f}", flush=True)
    print(f"\n{'method':>22} {'n_det':>6} {'CP bias':>10} {'+/-':>8}", flush=True)
    rows = []  # (method, n_det, bias, err)

    def record(m, nd, b, e):
        rows.append((m, nd, b, e))
        print(f"{m:>22} {nd:>6} {b:>+10.4f} {e:>8.4f}", flush=True)

    # --- free determinant ---
    ef, errf = run_single(a.lx, a.ly, a.nup, a.ndn, a.t, a.U, a.dt, a.nw,
                          a.nequil, a.nmeas, a.seed, "fixed")
    record("free-det", 1, ef - E0, errf)

    # --- adaptive (single self-consistent determinant) ---
    ea, erra = run_single(a.lx, a.ly, a.nup, a.ndn, a.t, a.U, a.dt, a.nw,
                          a.nequil, a.nmeas, a.seed, "adaptive")
    record("adaptive", 1, ea - E0, erra)

    # --- CASCI at several natural-orbital active windows (uses same W) ---
    cas_windows = []
    for nc in range(0, min(a.nup, a.ndn) + 1):
        for na in range(1, n - nc + 1):
            if a.nup - nc <= na and a.ndn - nc <= na:    # active electrons fit
                cas_windows.append((nc, na))
    cas_ncfgs = []
    for (nc, na) in sorted(set(cas_windows), key=lambda x: (x[0], x[1])):
        try:
            E0cas, ups, dns, coef, ncfg = casci_trial(K, a.U, W, a.nup, a.ndn, nc, na)
        except Exception:
            continue
        if ncfg in cas_ncfgs:
            continue
        cas_ncfgs.append(ncfg)
        e, err = run_multidet(a.lx, a.ly, a.nup, a.ndn, a.t, a.U, a.dt, a.nw,
                              a.nequil, a.nmeas, a.seed, ups, dns, coef)
        record(f"CASCI(core{nc},act{na})", ncfg, e - E0, err)

    # --- evolutionary selected-CI: one evolution, truncate to matched n_det ---
    k_max = min(len(_full(n, a.nup, a.ndn)), max(cas_ncfgs + [16]))
    rng = np.random.default_rng(a.seed)
    _, arch, civ = evolve_trial(h1, Up, n, a.nup, a.ndn, k_target=k_max, rng=rng,
                                ngen=300, n_candidates=96)
    order = np.argsort(-np.abs(civ)); arch = [arch[i] for i in order]; civ = civ[order]
    ga_ks = sorted(set([1] + cas_ncfgs))
    for k in ga_ks:
        if k > len(arch):
            continue
        ups, dns, coef = evolved_multidet(W, arch[:k], civ[:k])
        e, err = run_multidet(a.lx, a.ly, a.nup, a.ndn, a.t, a.U, a.dt, a.nw,
                              a.nequil, a.nmeas, a.seed, ups, dns, coef)
        record(f"evolution(k={k})", k, e - E0, err)



def _full(n, nup, ndn):
    from itertools import combinations
    return [(u, d) for u in combinations(range(n), nup)
                    for d in combinations(range(n), ndn)]


if __name__ == "__main__":
    main()
