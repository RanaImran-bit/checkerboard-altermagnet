#!/usr/bin/env python3
"""Phase 6 WS3 (a): VALIDATE the multi-determinant CPMC engine with a known-good
trial. Build the trial from the top-N largest-amplitude configurations of the
EXACT ED ground state (a truncated CI expansion), feed it to the verified md
machinery, and sweep N. This is a VALIDATION ONLY (needs ED -> small clusters);
it certifies the engine before the scalable natural-orbital CASCI (step b).

Self-validating: at N = all configurations the trial IS the exact ground state, so
the CPMC energy must equal E0 (zero constrained-path bias). As N grows from 1, the
bias must shrink monotonically toward 0 -- the definitive proof that a GOOD multi-
determinant trial reduces the CP bias (unlike the walker-sampled superposition).

    source tools/env.sh
    python pyqmc/validate_casci_ed.py --lx 2 --ly 2 --nup 2 --ndn 2 --U 4
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "ed"))


def det_from_config(sites, n, ne):
    """Slater-determinant orbital matrix (n x ne) in the SITE basis: columns are
    unit vectors at the occupied sites (sorted ascending)."""
    M = np.zeros((n, ne))
    for col, site in enumerate(sorted(sites)):
        M[site, col] = 1.0
    return M


def _hop_sign(occ, i, j):
    """Sign of <...|c^+_i c_j|...occ...> for a single-spin occupation set `occ`
    (sorted tuple): move a particle from occupied site j to empty site i. Sign =
    (-1)^(# occupied sites strictly between i and j) -- the Jordan-Wigner string,
    consistent with the sorted-ascending column convention of det_from_config."""
    lo, hi = (i, j) if i < j else (j, i)
    between = sum(1 for k in occ if lo < k < hi)
    return -1.0 if (between % 2) else 1.0


def ci_ground_state(K, pot_terms, n, nup, ndn):
    """Exact ground state in OUR sorted-column determinant convention: build the
    CI Hamiltonian by Slater-Condon (one-body K_ij hopping with JW sign; diagonal
    density-density interaction) over all configs, diagonalize. Returns E0 and a
    list of (coef, up_tuple, dn_tuple) sorted by |coef|. Self-consistent signs ->
    the full expansion IS the exact GS (full-N CPMC must give E0)."""
    from itertools import combinations
    ups = list(combinations(range(n), nup))
    dns = list(combinations(range(n), ndn))
    cfgs = [(u, d) for u in ups for d in dns]
    index = {c: m for m, c in enumerate(cfgs)}
    M = len(cfgs)
    H = np.zeros((M, M))
    Koff = [(i, j, K[i, j]) for i in range(n) for j in range(n)
            if i != j and abs(K[i, j]) > 1e-14]
    for m, (u, d) in enumerate(cfgs):
        us, ds = set(u), set(d)
        # diagonal: K_ii densities + density-density interaction
        diag = sum(K[i, i] for i in u) + sum(K[i, i] for i in d)
        for (a, sa, b, sb, V) in pot_terms:
            na = (a in us) if sa == 0 else (a in ds)
            nb = (b in us) if sb == 0 else (b in ds)
            diag += V * na * nb
        H[m, m] += diag
        # off-diagonal: one-body hopping K_ij c^+_i c_j (each spin)
        for (i, j, kij) in Koff:
            if j in us and i not in us:                      # up hop j->i
                nu = tuple(sorted((us - {j}) | {i}))
                mm = index[(nu, d)]
                H[mm, m] += kij * _hop_sign(u, i, j)
            if j in ds and i not in ds:                      # dn hop j->i
                nd = tuple(sorted((ds - {j}) | {i}))
                mm = index[(u, nd)]
                H[mm, m] += kij * _hop_sign(d, i, j)
    w, V = np.linalg.eigh(H)
    c0 = V[:, 0]
    out = sorted(((c0[m], cfgs[m][0], cfgs[m][1]) for m in range(M)),
                 key=lambda x: -abs(x[0]))
    return float(w[0]), out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=2); ap.add_argument("--ly", type=int, default=2)
    ap.add_argument("--nup", type=int, default=2); ap.add_argument("--ndn", type=int, default=2)
    ap.add_argument("--U", type=float, default=4.0); ap.add_argument("--t", type=float, default=1.0)
    ap.add_argument("--dt", type=float, default=0.01); ap.add_argument("--nw", type=int, default=300)
    ap.add_argument("--nequil", type=int, default=150); ap.add_argument("--nmeas", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--Ns", type=int, nargs="+", default=None,
                    help="list of #determinants to test (default: 1,2,4,8,...,all)")
    a = ap.parse_args()
    n = a.lx * a.ly

    from cpqmc import CPMC, square_hopping
    K = square_hopping(a.lx, a.ly, a.t)
    pot = [(i, 0, i, 1, a.U) for i in range(n)] if a.U != 0 else []
    E0, configs = ci_ground_state(K, pot, n, a.nup, a.ndn)
    ncfg = len(configs)
    # cross-check our CI E0 against QuSpin ED
    try:
        from hubbard_ed import build
        _, H, _, _ = build(a.lx, a.ly, a.nup, a.ndn, a.t, a.U)
        E0q = float(np.linalg.eigvalsh(H.toarray())[0])
        print(f"# CI E0 = {E0:.6f}  (QuSpin ED {E0q:.6f}, diff {abs(E0-E0q):.2e})")
    except Exception:
        pass
    print(f"# single-band {a.lx}x{a.ly}, nup={a.nup} ndn={a.ndn}, U={a.U}; "
          f"E0 = {E0:.6f}, {ncfg} configs")
    Ns = a.Ns or sorted(set([1, 2, 4, 8, 16, 32, 64, ncfg]) | {ncfg})
    Ns = [N for N in Ns if N <= ncfg]
    print(f"\n{'N_det':>6} {'CPMC E':>11} {'+/-':>8} {'dE(ED)':>9} {'trial-capture |<T|GS>|^2':>24}")
    for N in Ns:
        sub = configs[:N]
        coef = np.array([c for c, _, _ in sub])
        ups = [det_from_config(u, n, a.nup) for _, u, _ in sub]
        dns = [det_from_config(d, n, a.ndn) for _, _, d in sub]
        capture = float((coef ** 2).sum())          # fraction of GS captured
        coef = coef / np.linalg.norm(coef)
        q = CPMC(a.lx, a.ly, a.nup, a.ndn, t=a.t, U=a.U, dt=a.dt, nwalkers=a.nw,
                 seed=a.seed, trial="multidet")
        q.trial.set_multidet(ups, dns, coef)
        q.trial._refresh_olp(q.walkers)
        if q.walkers.olp[0] < 0:                      # orient trial sign to walkers
            q.trial.coef = -q.trial.coef
            q.trial._refresh_olp(q.walkers)
        e, err = q.run(nequil=a.nequil, nmeas=a.nmeas)
        print(f"{N:>6} {e:>11.4f} {err:>8.4f} {e - E0:>+9.4f} {capture:>24.5f}")


if __name__ == "__main__":
    main()
