#!/usr/bin/env python3
"""Algorithm A(b) — self-consistent EVOLUTIONARY trial with WALKER INJECTION.

This is the literal realization of "walkers as parents" (design note 00, idea 1):
the live CPMC walker population is mined for its dominant configurations, which are
injected as fresh genes into the evolutionary archive. The walkers explore the
correlated state under imaginary-time propagation, so they surface configurations
(e.g. pairing excitations) that local crossover/mutation from a poor reference
might never reach. The loop is:

    seed archive -> evolve (A(a)) -> install trial -> run CPMC window ->
    extract walkers' dominant configs -> inject as genes -> re-evolve -> repeat.

Plus a PAIRING SEED hook (idea for the d-wave target): plant paired double
excitations in the initial gene pool so the pairing manifold is reachable from
gen 0, independent of what the walkers find.

Built on evo_trial.py (A(a)); reuses the same selected-CI engine + ED gates. The
two-population firewall still holds: estimators come from the CPMC walkers; the
archive only defines the trial/constraint.
"""
from __future__ import annotations
import os, sys
from collections import Counter
from itertools import combinations
import numpy as np

_HERE = os.path.dirname(__file__)
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, "..", "..", "pyqmc"))
sys.path.insert(0, os.path.join(_HERE, "..", "..", "ed"))

from evo_trial import evolve_trial, evolved_multidet           # noqa: E402
from casci import transform_hubbard                            # noqa: E402
from cpqmc import CPMC, square_hopping                         # noqa: E402


# --------------------------------------------------------------------------
# Walker -> config extraction (the "parent" genes)
# --------------------------------------------------------------------------
def _walker_config(phi, W, ne):
    """Dominant configuration of a single-spin Slater determinant phi (n x ne) in
    the orbital basis W: the ne orbitals with the largest occupation in the
    walker's 1-RDM rho_orb = (W^T phi)(phi^T phi)^{-1}(W^T phi)^T."""
    A = W.T @ phi                                   # n_orb x ne
    rho = A @ np.linalg.solve(phi.T @ phi, A.T)
    occ = np.diag(rho).real
    return tuple(sorted(int(o) for o in np.argsort(-occ)[:ne]))


def extract_walker_configs(cpmc, W, nup, ndn, max_configs=48):
    """Weighted dominant configs across the walker ensemble. Returns up to
    `max_configs` (up_tuple, dn_tuple) ordered by total walker weight -- the
    configurations the population is actually visiting."""
    wk = cpmc.walkers
    cnt = Counter()
    for i in range(wk.nw):
        if wk.w[i] <= 0:
            continue
        cu = _walker_config(wk.phi_up[i], W, nup)
        cd = _walker_config(wk.phi_dn[i], W, ndn)
        cnt[(cu, cd)] += float(wk.w[i])
    return [c for c, _ in cnt.most_common(max_configs)]


# --------------------------------------------------------------------------
# Pairing seed (d-wave target): paired double excitations from the reference
# --------------------------------------------------------------------------
def pairing_seed_configs(n, nup, ndn, n_pairs=12):
    """Paired double excitations: simultaneously promote an up electron o->v and a
    down electron o->v (same orbital pair), the occupation-basis signature of a
    Cooper-pair excitation. Basis-agnostic (no momentum labels needed); the GA's
    energy fitness then keeps whichever paired excitations the model favors."""
    ref_u = list(range(nup)); ref_d = list(range(ndn))
    occ = sorted(set(ref_u) & set(ref_d))           # doubly-occupied orbitals
    virt = [o for o in range(n) if o not in set(ref_u) | set(ref_d)]
    out = []
    for o in reversed(occ):                          # excite from the top of the sea
        for v in virt:
            cu = tuple(sorted(set(ref_u) - {o} | {v}))
            cd = tuple(sorted(set(ref_d) - {o} | {v}))
            if len(cu) == nup and len(cd) == ndn:
                out.append((cu, cd))
            if len(out) >= n_pairs:
                return out
    return out


# --------------------------------------------------------------------------
# The self-consistent loop
# --------------------------------------------------------------------------
def _run_window(lx, ly, nup, ndn, t, U, dt, nw, seed, W, arch, civ,
                nequil, nmeas, K_dn=None, uxy=0.0, v=0.0):
    """Install the evolved trial, re-init walkers on its leading determinant, run a
    CPMC window. Returns (energy, err, cpmc)."""
    ups, dns, coef = evolved_multidet(W, arch, civ)
    q = CPMC(lx, ly, nup, ndn, t=t, U=U, dt=dt, nwalkers=nw, seed=seed,
             trial="multidet", K_dn=K_dn, uxy=uxy, v=v)
    q.trial.set_multidet(ups, dns, coef)
    lu = np.asarray(ups[0], float); ld = np.asarray(dns[0], float)
    q.walkers.phi_up = np.stack([lu.copy() for _ in range(q.walkers.nw)])
    q.walkers.phi_dn = np.stack([ld.copy() for _ in range(q.walkers.nw)])
    q.walkers.w = np.ones(q.walkers.nw)
    q.trial._refresh_olp(q.walkers)
    if q.walkers.olp[0] < 0:
        q.trial.coef = -q.trial.coef; q.trial._refresh_olp(q.walkers)
    e, err = q.run(nequil=nequil, nmeas=nmeas)
    return e, err, q


def scf_evolve(lx, ly, nup, ndn, t=1.0, U=4.0, dt=0.01, nw=400, seed=1,
               k_target=16, n_outer=4, seed_pairing=False,
               nequil=120, nmeas=200, K=None, K_dn=None, verbose=True):
    """Self-consistent evolutionary trial with walker injection. Returns
    (E_final, err_final, archive, civec, history)."""
    n = lx * ly
    Kmat = square_hopping(lx, ly, t) if K is None else np.asarray(K, float)
    _, W = np.linalg.eigh(Kmat)                      # K-eigenbasis (natural-like)
    h1, Up = transform_hubbard(Kmat, U, W)
    rng = np.random.default_rng(seed)

    seeds = [(tuple(range(nup)), tuple(range(ndn)))]
    if seed_pairing:
        seeds += pairing_seed_configs(n, nup, ndn)
    Eci, arch, civ = evolve_trial(h1, Up, n, nup, ndn, k_target=k_target, rng=rng,
                                  ngen=200, n_candidates=96, seed_configs=seeds)
    history = []
    e = err = float("nan")
    for outer in range(n_outer):
        e, err, q = _run_window(lx, ly, nup, ndn, t, U, dt, nw, seed + outer, W,
                                arch, civ, nequil, nmeas, K_dn=K_dn)
        wconfigs = extract_walker_configs(q, W, nup, ndn)
        # inject walker genes + keep current archive, re-evolve
        pool = list({*arch, *wconfigs})
        Eci, arch, civ = evolve_trial(h1, Up, n, nup, ndn, k_target=k_target,
                                      rng=rng, ngen=200, n_candidates=96,
                                      seed_configs=pool)
        n_new = len(set(wconfigs) - set(pool[:len(arch)]))
        history.append((outer, e, err, Eci, len(wconfigs)))
        if verbose:
            print(f"  outer {outer}: CPMC E={e:.4f}+/-{err:.4f}  "
                  f"trial E_CI={Eci:.4f}  injected {len(wconfigs)} walker configs",
                  flush=True)
    # final measurement with the converged trial
    e, err, _ = _run_window(lx, ly, nup, ndn, t, U, dt, nw, seed + n_outer, W,
                            arch, civ, nequil, nmeas, K_dn=K_dn)
    return e, err, arch, civ, history


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=2); ap.add_argument("--ly", type=int, default=2)
    ap.add_argument("--nup", type=int, default=2); ap.add_argument("--ndn", type=int, default=2)
    ap.add_argument("--U", type=float, default=4.0); ap.add_argument("--t", type=float, default=1.0)
    ap.add_argument("--k", type=int, default=16); ap.add_argument("--outer", type=int, default=4)
    ap.add_argument("--nw", type=int, default=400); ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--pairing", action="store_true")
    a = ap.parse_args()
    e, err, arch, civ, hist = scf_evolve(a.lx, a.ly, a.nup, a.ndn, t=a.t, U=a.U,
                                         nw=a.nw, seed=a.seed, k_target=a.k,
                                         n_outer=a.outer, seed_pairing=a.pairing)
    print(f"\nfinal: E={e:.5f} +/- {err:.5f}  |archive|={len(arch)}", flush=True)
