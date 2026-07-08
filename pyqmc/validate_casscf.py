#!/usr/bin/env python3
"""Phase 6 step (b): self-consistent natural-orbital CASCI trial (mini-CASSCF),
NO ED. The scalable multi-determinant trial the multidet engine needs.

Loop:
  1. run CPMC with the current trial (start: fixed free-electron GS);
  2. form the walker-ensemble charge 1-RDM, diagonalize -> natural orbitals +
     occupations (trial-agnostic);
  3. partition: freeze the top n_core orbitals, correlate the next n_active (CAS),
     drop the rest; transform the Hubbard integrals into that basis and diagonalize
     the small CAS Hamiltonian (Slater-Condon, our determinant convention);
  4. install the CASCI expansion as a frozen multi-determinant trial;
  5. iterate to self-consistency; final long CPMC run for the energy.

Acceptance gates (vs QuSpin ED):
  * FULL active space (n_core=0, n_active=n) MUST reproduce E0 (it is exact CI in
    a rotated basis) -> certifies the active-space + transform + CI path.
  * a REDUCED active space MUST reduce the constrained-path bias vs the fixed
    free-electron trial (the point of the method).

    source tools/env.sh
    python pyqmc/validate_casscf.py --lx 2 --ly 2 --nup 2 --ndn 2 --U 4 --ncore 1 --nact 2
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "ed"))
from cpqmc import CPMC, square_hopping                       # noqa: E402
from casci import (ensemble_rdm, natural_orbitals, casci_trial,  # noqa: E402
                   auto_active_space, occupation_entropy, cas_accuracy_estimate)


def ed_energy(lx, ly, nup, ndn, t, U):
    """Exact GS energy. Prefer QuSpin (sparse, scales); fall back to our own
    full-CI (casci_solve full active) when QuSpin is unavailable and the cluster
    is small enough. Pass --ed to skip entirely (e.g. on machines without QuSpin)."""
    try:
        from hubbard_ed import build
        _, H, _, _ = build(lx, ly, nup, ndn, t, U)
        return float(np.linalg.eigvalsh(H.toarray())[0])
    except Exception:
        from casci import casci_solve, transform_hubbard
        n = lx * ly
        h1, Up = transform_hubbard(square_hopping(lx, ly, t), U, np.eye(n))
        E0, _ = casci_solve(h1, Up, n, nup, ndn)
        return E0


def run_fixed(lx, ly, nup, ndn, t, U, nw, nequil, nmeas, seed):
    q = CPMC(lx, ly, nup, ndn, t=t, U=U, dt=0.01, nwalkers=nw, seed=seed, trial="fixed")
    return q.run(nequil=nequil, nmeas=nmeas)


def scf_casci(lx, ly, nup, ndn, t, U, n_core, n_active, nw, seed,
              n_iter=3, nequil=120, nmeas=200, verbose=True):
    """Self-consistent CASCI trial. Returns (E, err, E0_cas_final)."""
    n = lx * ly
    K = square_hopping(lx, ly, t)
    # iteration 0: a fixed-trial CPMC to seed the natural orbitals
    q = CPMC(lx, ly, nup, ndn, t=t, U=U, dt=0.01, nwalkers=nw, seed=seed, trial="fixed")
    q.run(nequil=nequil, nmeas=0)
    E = err = E0cas = float("nan")
    for it in range(n_iter):
        ru, rd = ensemble_rdm(q)
        occ, W = natural_orbitals(ru + rd)                   # charge natural orbitals
        E0cas, ups, dns, coef, ncfg = casci_trial(K, U, W, nup, ndn, n_core, n_active)
        # fresh ensemble on the new trial; install + orient sign to the walkers
        q = CPMC(lx, ly, nup, ndn, t=t, U=U, dt=0.01, nwalkers=nw, seed=seed + it + 1,
                 trial="multidet")
        q.trial.set_multidet(ups, dns, coef); q.trial._refresh_olp(q.walkers)
        if q.walkers.olp[0] < 0:
            q.trial.coef = -q.trial.coef; q.trial._refresh_olp(q.walkers)
        E, err = q.run(nequil=nequil, nmeas=nmeas)
        if verbose:
            print(f"  iter {it}: CAS({n_active}e? in {n_active}o, core {n_core}) "
                  f"E0_cas={E0cas:+.4f}  CPMC E={E:.4f} +/- {err:.4f}  "
                  f"occ={np.round(occ, 2)}")
    return E, err, E0cas


def scf_casci_auto(lx, ly, nup, ndn, t, U, nw, seed, n_iter=3,
                   nequil=120, nmeas=200, thr=0.05, max_active=8, verbose=True):
    """Self-consistent CASCI with an AUTOMATIC active space: each iteration the
    active orbitals are chosen from the natural-orbital occupations (fractional
    middle), so the CAS grows where the wavefunction is multireference. Reports an
    ED-free accuracy estimate (CAS truncation error) and the occupation entropy.
    Returns (E, err, n_active_final, dE_trunc, entropy)."""
    n = lx * ly
    K = square_hopping(lx, ly, t)
    q = CPMC(lx, ly, nup, ndn, t=t, U=U, dt=0.01, nwalkers=nw, seed=seed, trial="fixed")
    q.run(nequil=nequil, nmeas=0)
    E = err = dE = S = float("nan"); n_active = 0
    for it in range(n_iter):
        ru, rd = ensemble_rdm(q)
        occ, W = natural_orbitals(ru + rd)
        S = occupation_entropy(occ)
        n_core, n_active = auto_active_space(occ, nup, ndn, thr=thr, max_active=max_active)
        E0cas, dE, capture = cas_accuracy_estimate(K, U, W, nup, ndn, n_core, n_active)
        _, ups, dns, coef, ncfg = casci_trial(K, U, W, nup, ndn, n_core, n_active)
        q = CPMC(lx, ly, nup, ndn, t=t, U=U, dt=0.01, nwalkers=nw, seed=seed + it + 1,
                 trial="multidet")
        q.trial.set_multidet(ups, dns, coef); q.trial._refresh_olp(q.walkers)
        if q.walkers.olp[0] < 0:
            q.trial.coef = -q.trial.coef; q.trial._refresh_olp(q.walkers)
        E, err = q.run(nequil=nequil, nmeas=nmeas)
        if verbose:
            print(f"  iter {it}: entropy={S:.3f} -> CAS(core {n_core}, active {n_active}, "
                  f"{ncfg} dets)  E={E:.4f}+/-{err:.4f}  accuracy~|dE_trunc|={dE:.4f}")
    return E, err, n_active, dE, S


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=2); ap.add_argument("--ly", type=int, default=2)
    ap.add_argument("--nup", type=int, default=2); ap.add_argument("--ndn", type=int, default=2)
    ap.add_argument("--U", type=float, default=4.0); ap.add_argument("--t", type=float, default=1.0)
    ap.add_argument("--ncore", type=int, default=0); ap.add_argument("--nact", type=int, default=None)
    ap.add_argument("--nw", type=int, default=300); ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--iters", type=int, default=3)
    ap.add_argument("--nequil", type=int, default=120); ap.add_argument("--nmeas", type=int, default=200)
    ap.add_argument("--ed", type=float, default=None, help="pass exact E0 (skip QuSpin/CI)")
    ap.add_argument("--auto", action="store_true", help="auto-select n_active from occupations")
    ap.add_argument("--occthr", type=float, default=0.05, help="fractional-occupation threshold")
    ap.add_argument("--maxact", type=int, default=8, help="cap on auto n_active (cost)")
    a = ap.parse_args()
    n = a.lx * a.ly
    nact = a.nact if a.nact is not None else n - a.ncore

    E_ed = a.ed if a.ed is not None else ed_energy(a.lx, a.ly, a.nup, a.ndn, a.t, a.U)
    print(f"# {a.lx}x{a.ly} nup={a.nup} ndn={a.ndn} U={a.U}  ED E0 = {E_ed:.4f}")

    ef, errf = run_fixed(a.lx, a.ly, a.nup, a.ndn, a.t, a.U, a.nw, a.nequil, a.nmeas, a.seed)
    print(f"FIXED free-electron trial : E={ef:.4f} +/- {errf:.4f}  dE(ED)={ef - E_ed:+.4f}")

    if a.auto:
        print(f"CASCI-SCF AUTO (occthr={a.occthr}, maxact={a.maxact}):")
        ec, errc, nact_f, dE, S = scf_casci_auto(a.lx, a.ly, a.nup, a.ndn, a.t, a.U,
                                                 a.nw, a.seed, a.iters, a.nequil, a.nmeas,
                                                 a.occthr, a.maxact)
        print(f"CASCI-SCF AUTO trial     : E={ec:.4f} +/- {errc:.4f}  dE(ED)={ec - E_ed:+.4f}  "
              f"n_active={nact_f}  entropy={S:.3f}  accuracy~|dE_trunc|={dE:.4f}")
        print(f"  fixed |dE|={abs(ef - E_ed):.4f} -> AUTO-CASCI |dE|={abs(ec - E_ed):.4f}  "
              f"(predicted ceiling ~{dE:.3f})")
        return

    print(f"CASCI-SCF (n_core={a.ncore}, n_active={nact}):")
    ec, errc, e0 = scf_casci(a.lx, a.ly, a.nup, a.ndn, a.t, a.U, a.ncore, nact,
                             a.nw, a.seed, a.iters, a.nequil, a.nmeas)
    print(f"CASCI-SCF trial          : E={ec:.4f} +/- {errc:.4f}  dE(ED)={ec - E_ed:+.4f}  "
          f"(CAS CI E0={e0:.4f})")
    tag = "FULL active -> must equal ED" if nact == n - a.ncore and a.ncore == 0 else "reduced active"
    print(f"  [{tag}]  fixed |dE|={abs(ef - E_ed):.4f} -> CASCI |dE|={abs(ec - E_ed):.4f}")


if __name__ == "__main__":
    main()
