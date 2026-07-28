#!/usr/bin/env python3
"""QuSpin-FREE ED gate for the checkerboard dynamic pairing susceptibility VERTEX
(s, dx2-y2, dxy). numpy + scipy only -> runs on qmc48 (no QuSpin). Same purpose as
validate_checkerboard_chid.py: confirm the SIGN of the connected-vertex chi vs the
CP-AFQMC, so we know whether to trust the dxy vertex.

Method: sector product basis (up_cfgs x dn_cfgs) from checkerboard_ed_np. The pair
operator Delta^dag = sum f_a(delta) c^+_{m up} c^+_{m+delta dn} and the single-particle
c^+_a are built as rectangular blocks with WITHIN-spin Jordan-Wigner signs; the
cross-spin (-1)^nup sign cancels because C = |<n|Delta^dag|0>|^2 is quadratic. Full-Fock
dense diag => keep <= 6 sites (2x3 / 3x2).

  python3 pyqmc/validate_checkerboard_chid_np.py --lx 3 --ly 2 --nup 2 --ndn 2 --U 4 --delta 0.2
"""
from __future__ import annotations
import os, sys, argparse
import numpy as np
from scipy.sparse import csr_matrix, identity, kron, diags
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "Checkerboard_Model"))
from cpqmc import CPMC
import checkerboard as cb
from checkerboard_ed import build_hopping
from checkerboard_ed_np import spin_configs, hop_block, interaction_diag, _popcount_below


def create_block(cfgs_from, cfgs_to, m):
    """Rectangular c^+_m : sector(np) -> sector(np+1), within-spin JW sign."""
    rows, cols, vals = [], [], []
    for c in range(len(cfgs_from)):
        mm = int(cfgs_from[c])
        if (mm >> m) & 1:                         # m must be empty
            continue
        s = _popcount_below(mm, m)
        r = int(np.searchsorted(cfgs_to, np.uint64(mm | (1 << m))))
        rows.append(r); cols.append(c); vals.append(1.0 if (s & 1) == 0 else -1.0)
    return csr_matrix((vals, (rows, cols)), shape=(len(cfgs_to), len(cfgs_from)))


def sector_H(N, terms, up, dn, U):
    Hu = hop_block(up, terms); Hd = hop_block(dn, terms)
    D = (interaction_diag(up, dn, U).reshape(len(up), len(dn)) if U != 0
         else np.zeros((len(up), len(dn))))
    return (kron(Hu, identity(len(dn)), format="csr")
            + kron(identity(len(up)), Hd, format="csr")
            + diags(D.ravel(), 0, format="csr")).toarray()


def ed_cb_chid_np(lx, ly, nup, ndn, t0, t1, t2, U, taus, Ffac):
    N = lx * ly
    tk = build_hopping(lx, ly, t0, t1, t2)
    terms = [(i, j, tk[i, j]) for i in range(N) for j in range(N) if i != j and abs(tk[i, j]) > 1e-15]
    up0 = spin_configs(N, nup);   dn0 = spin_configs(N, ndn)
    up1 = spin_configs(N, nup + 1); dn1 = spin_configs(N, ndn + 1)
    Cu0, Cd0 = len(up0), len(dn0)

    # ground state in (nup,ndn)
    H0 = sector_H(N, terms, up0, dn0, U)
    w0, v0 = np.linalg.eigh(H0); E0 = float(w0[0]); psi0 = v0[:, 0]

    # +pair sector (nup+1,ndn+1) -> full-susceptibility Lehmann
    H2 = sector_H(N, terms, up1, dn1, U)
    E2, V2 = np.linalg.eigh(H2)
    # +1up and +1dn sectors -> single-particle GFs (bubble)
    H1u = sector_H(N, terms, up1, dn0, U); E1u, V1u = np.linalg.eigh(H1u); dE1u = E1u - E0
    H1d = sector_H(N, terms, up0, dn1, U); E1d, V1d = np.linalg.eigh(H1d); dE1d = E1d - E0

    Aup = [create_block(up0, up1, m) for m in range(N)]   # c^+_{m up}: up0 -> up1
    Adn = [create_block(dn0, dn1, m) for m in range(N)]   # c^+_{m dn}: dn0 -> dn1
    Iu0 = identity(Cu0, format="csr"); Id0 = identity(Cd0, format="csr")

    # single-particle amplitudes M[a,n] = <n| c^+_{a} |0> in the +1 sectors
    Mu = np.zeros((N, len(E1u))); Md = np.zeros((N, len(E1d)))
    for a in range(N):
        cu = kron(Aup[a], Id0, format="csr").dot(psi0)     # in (nup+1,ndn) product basis
        cd = kron(Iu0, Adn[a], format="csr").dot(psi0)     # in (nup,ndn+1) product basis
        Mu[a, :] = V1u.T @ cu; Md[a, :] = V1d.T @ cd

    res = {}
    for tag, F in Ffac.items():
        # Delta^dag_a = sum_{m,j} F[m,j] c^+_{m up} c^+_{j dn}  (cross-spin sign cancels in |.|^2)
        Op = None
        for m in range(N):
            for j in range(N):
                f = float(F[m, j])
                if abs(f) < 1e-12:
                    continue
                term = kron(Aup[m], Adn[j], format="csr") * f
                Op = term if Op is None else (Op + term)
        vec = Op.dot(psi0)                                  # in (nup+1,ndn+1) product basis
        a = V2.T @ vec; dE = E2 - E0
        Cfull = np.array([float(np.sum(a ** 2 * np.exp(-d * dE))) for d in taus])
        Cbub = np.empty(len(taus))
        for k, tau in enumerate(taus):
            Gu = (Mu * np.exp(-dE1u * tau)) @ Mu.T
            Gd = (Md * np.exp(-dE1d * tau)) @ Md.T
            Cbub[k] = float(np.sum(Gu * (F @ Gd @ F.T)))
        res[tag] = (Cfull, Cbub, Cfull - Cbub)
    return E0, res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=3); ap.add_argument("--ly", type=int, default=2)
    ap.add_argument("--nup", type=int, default=2); ap.add_argument("--ndn", type=int, default=2)
    ap.add_argument("--U", type=float, default=4.0)
    ap.add_argument("--t0", type=float, default=-1.0); ap.add_argument("--t1", type=float, default=0.3)
    ap.add_argument("--delta", type=float, default=0.2)
    ap.add_argument("--dt", type=float, default=0.05); ap.add_argument("--nw", type=int, default=240)
    ap.add_argument("--nequil", type=int, default=60); ap.add_argument("--nblocks", type=int, default=80)
    ap.add_argument("--bp", type=int, default=20)
    a = ap.parse_args()
    if a.lx * a.ly > 6:
        print(f"WARNING: {a.lx*a.ly} sites; dense sector diag gets big beyond 6.")
    t2 = -a.delta
    K = cb.checkerboard_hopping(a.lx, a.ly, a.t0, a.t1, t2)
    Fs, Fd = cb.nn_bond_factors(a.lx, a.ly); Fdxy = cb.diag_bond_factors(a.lx, a.ly)
    Ffac = {"s": Fs, "d": Fd, "dxy": Fdxy}

    q = CPMC(a.lx, a.ly, a.nup, a.ndn, U=a.U, dt=a.dt, nwalkers=a.nw, seed=1, K=K, K_dn=None)
    qm = cb.run_bp_chid_cb(q, Ffac, nequil=a.nequil, nblocks=a.nblocks, bp=a.bp)
    taus = np.array(qm["taus"])
    E0, ed = ed_cb_chid_np(a.lx, a.ly, a.nup, a.ndn, a.t0, a.t1, t2, a.U, taus, Ffac)

    def winint(C): return a.dt * (C[1:-1].sum() + 0.5 * (C[0] + C[-1]))
    print(f"# checkerboard chi ED gate (QuSpin-free)  L={a.lx}x{a.ly}  n={2*a.nup/(a.lx*a.ly):.3f}  U={a.U}  delta={a.delta}")
    print(f"# E0(ED)={E0:.5f}   window 0..{taus[-1]:.2f}")
    print(f"\n### EQUAL-TIME vertex C^vtx(tau=0) -- the quantity Figs 6-8 sum ###")
    print(f"{'ch':5} {'ED':>11} {'QMC':>11}  sign")
    et_ok = True
    for tag in ("s", "d", "dxy"):
        edE = ed[tag][2][0]                          # tau=0 slice of the ED vertex
        qE  = qm[f"C_{tag}_tau_vertex"][0]           # tau=0 slice of the CP-AFQMC vertex
        sgn = "OK" if (np.sign(edE) == np.sign(qE) or abs(edE) < 2e-3) else "FLIP"
        et_ok = et_ok and (sgn == "OK")
        print(f"{tag:5} {edE:>+11.4f} {qE:>+11.4f}  {sgn}")
    print(f"\n### tau-INTEGRATED susceptibility vertex chi^vtx -- the sound observable ###")
    print(f"{'ch':5} {'chiF ED':>10} {'chiF QMC':>10} | {'chiV ED':>11} {'chiV QMC':>11} {'+/-':>8}  sign")
    int_ok = True
    for tag in ("s", "d", "dxy"):
        Cf, Cb, Cv = ed[tag]
        edF = winint(Cf); edV = winint(Cv)
        qF = qm[f"chi_{tag}"]; qV = qm[f"chi_{tag}_vertex"]; qVe = qm[f"chi_{tag}_vertex_err"]
        sgn = "OK" if (np.sign(edV) == np.sign(qV) or abs(edV) < 2e-3) else "FLIP"
        int_ok = int_ok and (sgn == "OK")
        print(f"{tag:5} {edF:>10.4f} {qF:>10.4f} | {edV:>+11.4f} {qV:>+11.4f} {qVe:>8.4f}  {sgn}")
    print(f"\nVERDICT  equal-time vertex: {'sign-faithful' if et_ok else 'CP-BIASED (some sign FLIP)'}"
          f"  |  integrated susceptibility: {'sign-faithful' if int_ok else 'CP-BIASED'}")


if __name__ == "__main__":
    main()
