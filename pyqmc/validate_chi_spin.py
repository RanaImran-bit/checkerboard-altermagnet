#!/usr/bin/env python3
"""ED gate for the MOMENTUM-RESOLVED T=0 spin susceptibility (chi_zz(q) work item W1,
docs/PLAN_chi_spin_AHE.md): pyqmc CPMC run_bp_chi_spin vs exact Lehmann ED, single-band
(altermagnet) Hubbard, every q on the (lx,ly) grid.

    C_q(tau) = (1/N) sum_n |<n|S^z_q|0>|^2 e^{-(E_n-E_0) tau},
    S^z_q    = sum_i e^{-iq.r_i} S^z_i,  q = 2pi(kx/lx, ky/ly),  site i = x*ly + y,
    chi_q^win = trapezoid of C_q(tau) over the BP window [0, bp*dt]

(the same windowed one-sided T=0 Kubo convention as validate_chi; at (pi,pi) this
reduces to N * the staggered gate of validate_chi -- cross-checked here).

NB cross-engine: the finite-T engines (SpinDQMC / FTCPMC) use chi = int_0^beta dtau,
whose beta->inf limit is 2x this one-sided T=0 value: chi_finiteT -> 2 * chi_T0.
(Verified against spin_susc.ed_chi_spin at beta=12: equal-time C_q(0) identical,
chi ratio exactly 2.)

    source tools/env.sh
    python pyqmc/validate_chi_spin.py --lx 4 --ly 2 --nup 4 --ndn 4 --U 0     # exact gate
    python pyqmc/validate_chi_spin.py --lx 2 --ly 2 --nup 2 --ndn 2 --U 4
    python pyqmc/validate_chi_spin.py --lx 2 --ly 2 --nup 2 --ndn 2 --U 4 --tam 0.3
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMBA_NUM_THREADS", "1")
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
import argparse, json, sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))


def _sector_hop(K, n, npart):
    """Dense one-species hopping Hamiltonian H[a,b] = sum_ij K_ij <a|c^+_i c_j|b> on the
    fixed-particle-number bitmask basis (Jordan-Wigner sign = parity of occupied sites
    strictly between i and j). Returns (H, configs, occ) with occ[a,i] in {0,1}."""
    from itertools import combinations
    configs = [sum(1 << s for s in c) for c in combinations(range(n), npart)]
    index = {c: a for a, c in enumerate(configs)}
    D = len(configs)
    occ = np.array([[(c >> i) & 1 for i in range(n)] for c in configs], dtype=np.int8)
    H = np.zeros((D, D))
    for a, c in enumerate(configs):
        for i in range(n):
            for j in range(n):
                if abs(K[i, j]) < 1e-14:
                    continue
                if i == j:                              # diagonal (chemical/on-site) term
                    if (c >> i) & 1:
                        H[a, a] += K[i, i]
                    continue
                if not ((c >> j) & 1) or ((c >> i) & 1):
                    continue                            # need j occupied, i empty
                lo, hi = (i, j) if i < j else (j, i)
                between = c & (((1 << hi) - 1) ^ ((1 << (lo + 1)) - 1))
                sgn = -1.0 if bin(between).count("1") % 2 else 1.0
                H[index[(c ^ (1 << j)) | (1 << i)], a] += sgn * K[i, j]
    return H, configs, occ


def ed_chi_spin_T0(lx, ly, nup, ndn, t0, U, taus, tam=0.0, t1=0.0, tp=0.0):
    """T=0 Lehmann C_q(tau) on the full (lx,ly) q-grid. Numpy-only sector ED
    (fixed (nup,ndn); H = Hup x I + I x Hdn + U diag -- S^z_q is DIAGONAL in the
    occupation basis, so the Lehmann sum needs only the sector eigenbasis).
    Same site convention/i-ordering (i = x*ly + y) as the CPMC engine. Returns
    (Ctau_q[L,lx,ly], chi_q_full[lx,ly]) with chi_q_full the UNWINDOWED sum |a|^2/dE."""
    from cpqmc import am_hopping
    n = lx * ly
    Ku, Kd = am_hopping(lx, ly, t0, tam, t1, tp)
    Hu, cu, occu = _sector_hop(Ku, n, nup)
    Hd, cd, occd = _sector_hop(Kd, n, ndn)
    Du, Dd = Hu.shape[0], Hd.shape[0]
    H = (np.kron(Hu, np.eye(Dd)) + np.kron(np.eye(Du), Hd)
         + np.diag(U * (occu[:, None, :] * occd[None, :, :]).sum(-1).ravel().astype(float)))
    E, V = np.linalg.eigh(H)
    psi0 = V[:, 0]; dE = E - E[0]
    xm = np.arange(n) // ly; ym = np.arange(n) % ly
    mz = 0.5 * (occu[:, None, :] - occd[None, :, :]).reshape(Du * Dd, n)  # <a,b|S^z_i|a,b>
    L = len(taus)
    Ctau_q = np.zeros((L, lx, ly)); chi_q = np.zeros((lx, ly))
    for kx in range(lx):
        for ky in range(ly):
            ph = np.exp(-2j * np.pi * (kx * xm / lx + ky * ym / ly))
            o = mz @ ph                                  # diagonal of S^z_q
            a2 = np.abs(V.T @ (o * psi0)) ** 2           # |<m|S^z_q|0>|^2 (V real)
            Ctau_q[:, kx, ky] = [float(np.sum(a2 * np.exp(-d * dE))) for d in taus]
            mask = dE > 1e-9
            chi_q[kx, ky] = float(np.sum(a2[mask] / dE[mask]))
    return Ctau_q / n, chi_q / n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=2); ap.add_argument("--ly", type=int, default=2)
    ap.add_argument("--nup", type=int, default=2); ap.add_argument("--ndn", type=int, default=2)
    ap.add_argument("--U", type=float, default=4.0); ap.add_argument("--t", type=float, default=1.0)
    ap.add_argument("--tam", type=float, default=0.0); ap.add_argument("--t1", type=float, default=0.0)
    ap.add_argument("--tp", type=float, default=0.0)
    ap.add_argument("--dt", type=float, default=0.01); ap.add_argument("--nw", type=int, default=300)
    ap.add_argument("--nequil", type=int, default=200); ap.add_argument("--nblocks", type=int, default=40)
    ap.add_argument("--bp", type=int, default=20); ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--tol-nsig", type=float, default=3.0, help="pass if |dev| < max(nsig*err, rtol*ED)")
    ap.add_argument("--rtol", type=float, default=0.05)
    ap.add_argument("-o", "--out", help="write JSON record")
    a = ap.parse_args()

    from cpqmc import CPMC, am_hopping
    Ku, Kd = am_hopping(a.lx, a.ly, a.t, a.tam, a.t1, a.tp)
    q = CPMC(a.lx, a.ly, a.nup, a.ndn, t=a.t, U=a.U, dt=a.dt, nwalkers=a.nw,
             seed=a.seed, K=Ku, K_dn=Kd)
    qm = q.run_bp_chi_spin(nequil=a.nequil, nblocks=a.nblocks, bp=a.bp)
    taus = np.array(qm["taus"])
    Cq = np.array(qm["Ctau_q"]); Cqe = np.array(qm["Cerr_q"])          # (L, lx, ly)
    chi_q = np.array(qm["chi_q"]); chi_qe = np.array(qm["chi_q_err"])

    Ce, chi_ed_full = ed_chi_spin_T0(a.lx, a.ly, a.nup, a.ndn, a.t, a.U, taus,
                                     a.tam, a.t1, a.tp)
    # windowed ED chi (same trapezoid as the QMC)
    chi_ed_win = a.dt * (Ce[1:-1].sum(axis=0) + 0.5 * (Ce[0] + Ce[-1]))

    print(f"# chi_zz(q) T=0 gate: {a.lx}x{a.ly} nup={a.nup} ndn={a.ndn} U={a.U} "
          f"tam={a.tam} t1={a.t1} tp={a.tp}  (bp={a.bp}, dt={a.dt}, window tau<={taus[-1]:.2f})")
    print(f"{'q=(kx,ky)':>10} {'ED chi_win':>11} {'CPMC chi':>10} {'+/-':>8} {'z':>6}  "
          f"{'ED C(0)':>9} {'CPMC C(0)':>10}")
    npass = ntot = 0
    for kx in range(a.lx):
        for ky in range(a.ly):
            ed, qc, qe = chi_ed_win[kx, ky], chi_q[kx, ky], chi_qe[kx, ky]
            dev = abs(qc - ed)
            tol = max(a.tol_nsig * qe, a.rtol * max(abs(ed), 1e-3))
            ok = dev < tol
            npass += ok; ntot += 1
            z = dev / qe if qe > 0 else float("nan")
            print(f"  ({kx},{ky})   {ed:>11.4f} {qc:>10.4f} {qe:>8.4f} {z:>6.2f}  "
                  f"{Ce[0, kx, ky]:>9.4f} {Cq[0, kx, ky]:>10.4f}"
                  + ("" if ok else "   << FAIL"))
    verdict = "PASS" if npass == ntot else "FAIL"
    print(f"\n{verdict}: {npass}/{ntot} q-points within max({a.tol_nsig} sigma, "
          f"{100*a.rtol:.0f}% ED)  [chi window 0..{taus[-1]:.2f}]")
    kpi = (a.lx // 2, a.ly // 2)
    print(f"(pi,pi) cross-check vs validate_chi convention: N*chi_q = "
          f"{a.lx*a.ly*chi_q[kpi]:.4f} (ED windowed {a.lx*a.ly*chi_ed_win[kpi]:.4f}, "
          f"ED full {a.lx*a.ly*chi_ed_full[kpi]:.4f})")
    if a.out:
        with open(a.out, "w") as f:
            json.dump({"params": vars(a), "qmc": qm,
                       "ed_chi_win": chi_ed_win.tolist(),
                       "ed_chi_full": chi_ed_full.tolist(),
                       "ed_Ctau_q": Ce.tolist(), "verdict": verdict}, f, indent=1)
        print(f"wrote {a.out}")
    sys.exit(0 if verdict == "PASS" else 1)


if __name__ == "__main__":
    main()
