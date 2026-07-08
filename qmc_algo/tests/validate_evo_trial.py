#!/usr/bin/env python3
"""ED-gated validation of Algorithm A(a) — the evolutionary selected-CI trial.

Gate 1 (EXACTNESS): with the archive allowed to grow to the full configuration
space, the evolved selected-CI energy must equal the exact ED ground state
(basis-invariant CI). NO physics claim survives if this fails.

Gate 2 (BIAS): feed the evolved archive (built with NO ED) as the CPMC trial and
sweep its size k. The constrained-path bias E_CPMC - E0 must fall monotonically
toward 0 as k grows -- the thesis of design note 00, and it must beat the single
free-electron determinant. We print the GA selected-CI energy alongside the
ED-optimal truncation (configs ranked by the EXACT CI coefficient) to show the
GA approaches the oracle ordering WITHOUT using ED.

    source ../../tools/env.sh   # (from this dir) or repo tools/env.sh
    python validate_evo_trial.py --lx 2 --ly 2 --nup 2 --ndn 2 --U 4
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np

_HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(_HERE, "..", "src"))
sys.path.insert(0, os.path.join(_HERE, "..", "..", "pyqmc"))
sys.path.insert(0, os.path.join(_HERE, "..", "..", "ed"))

from evo_trial import (evolve_trial, selected_ci_solve, evolved_multidet,  # noqa
                       _space_size)
from casci import transform_hubbard, build_multidet_trial
from cpqmc import CPMC, square_hopping
from validate_casci_ed import ci_ground_state


def ed_reference(lx, ly, nup, ndn, t, U):
    """Exact GS energy + ED-ranked configs (in the K-eigenbasis W) for the oracle
    truncation. Uses QuSpin if available, else the in-house CI."""
    K = square_hopping(lx, ly, t)
    n = lx * ly
    _, W = np.linalg.eigh(K)                     # natural-orbital-like basis
    h1, Upqrs = transform_hubbard(K, U, W)
    # ED-optimal expansion in basis W = full selected CI over all configs
    from itertools import combinations
    full = [(u, d) for u in combinations(range(n), nup)
                    for d in combinations(range(n), ndn)]
    E0, civec, _ = selected_ci_solve(full, h1, Upqrs)
    try:
        from hubbard_ed import build
        _, H, _, _ = build(lx, ly, nup, ndn, t, U)
        E0 = float(np.linalg.eigvalsh(H.toarray())[0])
    except Exception:
        pass
    order = np.argsort(-np.abs(civec))
    ranked = [full[i] for i in order]
    return K, W, h1, Upqrs, E0, full, civec, ranked


def run_cpmc_with_archive(lx, ly, nup, ndn, t, U, dt, nw, nequil, nmeas, seed,
                          W, archive, civec):
    ups, dns, coef = evolved_multidet(W, archive, civec)
    q = CPMC(lx, ly, nup, ndn, t=t, U=U, dt=dt, nwalkers=nw, seed=seed,
             trial="multidet")
    q.trial.set_multidet(ups, dns, coef)
    # Re-initialize the walkers from the trial's LEADING determinant. The trial
    # lives in the K-eigenbasis (same basis as the default free-electron walkers),
    # so those start orthogonal to it (overlap 0/1) -> overlap collapse / nan for
    # degenerate multideterminant trials. Seeding walkers on the leading
    # determinant guarantees a finite starting overlap; propagation then spreads.
    lead_up = np.asarray(ups[0], float); lead_dn = np.asarray(dns[0], float)
    q.walkers.phi_up = np.stack([lead_up.copy() for _ in range(q.walkers.nw)])
    q.walkers.phi_dn = np.stack([lead_dn.copy() for _ in range(q.walkers.nw)])
    q.walkers.w = np.ones(q.walkers.nw)
    q.trial._refresh_olp(q.walkers)
    if q.walkers.olp[0] < 0:                     # orient trial sign to walkers
        q.trial.coef = -q.trial.coef
        q.trial._refresh_olp(q.walkers)
    return q.run(nequil=nequil, nmeas=nmeas)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=2); ap.add_argument("--ly", type=int, default=2)
    ap.add_argument("--nup", type=int, default=2); ap.add_argument("--ndn", type=int, default=2)
    ap.add_argument("--U", type=float, default=4.0); ap.add_argument("--t", type=float, default=1.0)
    ap.add_argument("--dt", type=float, default=0.01); ap.add_argument("--nw", type=int, default=300)
    ap.add_argument("--nequil", type=int, default=150); ap.add_argument("--nmeas", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--ks", type=int, nargs="+", default=None)
    ap.add_argument("--no-cpmc", action="store_true", help="run only the exactness gate")
    a = ap.parse_args()
    n = a.lx * a.ly
    rng = np.random.default_rng(a.seed)

    K, W, h1, Upqrs, E0, full, civec_full, ranked = ed_reference(
        a.lx, a.ly, a.nup, a.ndn, a.t, a.U)
    ncfg = len(full)
    print(f"# single-band {a.lx}x{a.ly}, nup={a.nup} ndn={a.ndn}, U={a.U}, "
          f"t={a.t}; ED E0 = {E0:.6f}; {ncfg} configs total")

    # ---- Gate 1: exactness (archive -> full space) ----
    E_full, arch_full, _ = evolve_trial(h1, Upqrs, n, a.nup, a.ndn,
                                        k_target=ncfg, rng=rng, ngen=200,
                                        n_candidates=min(128, ncfg))
    dE = abs(E_full - E0)
    ok = dE < 1e-6
    print(f"\n[Gate 1 EXACTNESS] evolved full-space CI E0 = {E_full:.8f}  "
          f"|E - ED| = {dE:.2e}  -> {'PASS' if ok else 'FAIL'}  "
          f"(|archive|={len(arch_full)}/{ncfg})")
    if a.no_cpmc:
        return

    # ---- Gate 2: CP-bias sweep vs archive size k ----
    # Evolve ONCE to a generous size, then truncate by |CI coefficient| -> NESTED
    # truncations (so the CP bias is monotone), and the HF reference is free to be
    # dropped at small k (essential in a multireference state).
    ks = a.ks or sorted(set([1, 2, 4, 8, 16, min(32, ncfg), ncfg]))
    ks = [k for k in ks if k <= ncfg]
    k_max = min(ncfg, max(ks))
    Emax, arch_full_evo, civ_full_evo = evolve_trial(
        h1, Upqrs, n, a.nup, a.ndn, k_target=k_max, rng=rng, ngen=300,
        n_candidates=min(96, ncfg))
    order_evo = np.argsort(-np.abs(civ_full_evo))
    arch_evo = [arch_full_evo[i] for i in order_evo]      # descending |coef|
    civ_evo = civ_full_evo[order_evo]
    gs_amp = {c: civec_full[i] for i, c in enumerate(full)}   # exact GS amplitudes

    print(f"\n[Gate 2 BIAS] CPMC constrained-path bias vs evolved-trial size k "
          f"(single evolution to k_max={k_max}, E_CI={Emax:.4f})")
    print(f"{'k':>5} {'GA_CI_E':>11} {'oracle_E':>11} {'CPMC_E':>11} {'+/-':>8} "
          f"{'bias':>9} {'|<T|GS>|^2':>10}")
    # single free-electron determinant baseline (k=1 fixed trial)
    qf = CPMC(a.lx, a.ly, a.nup, a.ndn, t=a.t, U=a.U, dt=a.dt, nwalkers=a.nw,
              seed=a.seed, trial="fixed")
    ef, errf = qf.run(nequil=a.nequil, nmeas=a.nmeas)
    print(f"{'1*':>5} {'':>11} {'':>11} {ef:>11.4f} {errf:>8.4f} "
          f"{ef - E0:>+9.4f} {'(free-det)':>10}")
    for k in ks:
        arch = arch_evo[:k]
        civ = civ_evo[:k]
        Ek, _, _ = selected_ci_solve(arch, h1, Upqrs)            # variational E of subset
        Eo, _, _ = selected_ci_solve(ranked[:k], h1, Upqrs)      # ED-ranked oracle
        # true squared overlap of the (normalized) truncated trial with exact GS
        cn = civ / np.linalg.norm(civ)
        ov = sum(cn[i] * gs_amp.get(arch[i], 0.0) for i in range(k))
        e, err = run_cpmc_with_archive(a.lx, a.ly, a.nup, a.ndn, a.t, a.U, a.dt,
                                       a.nw, a.nequil, a.nmeas, a.seed, W, arch, civ)
        print(f"{k:>5} {Ek:>11.4f} {Eo:>11.4f} {e:>11.4f} {err:>8.4f} "
              f"{e - E0:>+9.4f} {ov * ov:>10.4f}")


if __name__ == "__main__":
    main()
