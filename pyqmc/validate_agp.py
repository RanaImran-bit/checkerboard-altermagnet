#!/usr/bin/env python3
"""ED gate for the AGP / number-projected-BCS geminal trial (pyqmc/agp.py).

Builds the EXACT many-body AGP state and a Slater walker in QuSpin on a tiny
cluster, then checks the closed-form overlap / Green's-function / local-energy
formulas in agp.py against the exact mixed expectation values:

    |Psi_AGP> = (sum_ij F_ij c^+_{i up} c^+_{j dn})^N |0>     (apply pair op N times)
    |Phi>     = Slater determinant with up/dn orbital matrices A, B

  overlap:  det(A^T F B)                  vs  <Psi_AGP|Phi>           (up to global sign)
  Gu[i,j]:  (F B M^{-1} A^T)_{ij}          vs  <Psi_AGP|c^+_{iu}c_{ju}|Phi>/<.|.>
  Gd[i,j]:  (F^T A M^{-T} B^T)_{ij}        vs  <Psi_AGP|c^+_{id}c_{jd}|Phi>/<.|.>
  E_local:  Wick(Gu,Gd) on H              vs  <Psi_AGP|H|Phi>/<Psi_AGP|Phi>

The Green's functions are pure ratios (normalization-independent), so a match to
~1e-10 unambiguously validates the formulas + the BdG->geminal construction. This
is the rigorous "generalized Slater/HFB overlap" the sister pph-qmc effort needs.

    source tools/env.sh
    python pyqmc/validate_agp.py --lx 2 --ly 2 --nup 1 --ndn 1 --U 4 --lam 0.15
    python pyqmc/validate_agp.py --lx 2 --ly 2 --nup 2 --ndn 2 --U 4 --t1 0.3 --lam 0.2
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMBA_NUM_THREADS", "1")
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
import argparse, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from cpqmc import am_hopping, LatticeModel
from agp import dwave_geminal, agp_contractions, agp_local_energy


def build_ops(lx, ly, Ku, Kd, U):
    """QuSpin full-Fock-space pieces: H (Hubbard with spin-dependent hopping Ku/Kd),
    per-site number diagonals, and a factory for creation/number operators."""
    from quspin.basis import spinful_fermion_basis_general
    from quspin.operators import hamiltonian
    n = lx * ly
    basis = spinful_fermion_basis_general(n)
    nc = dict(check_pcon=False, check_symm=False, check_herm=False)
    # H = sum_ij Ku_ij c^+_iu c_ju + Kd_ij c^+_id c_jd + U sum_i n_iu n_id
    up_hop = [[Ku[i, j], i, j] for i in range(n) for j in range(n) if abs(Ku[i, j]) > 1e-15]
    dn_hop = [[Kd[i, j], i, j] for i in range(n) for j in range(n) if abs(Kd[i, j]) > 1e-15]
    inter = [[U, i, i] for i in range(n)]
    static = [["+-|", up_hop], ["|+-", dn_hop], ["n|n", inter]]
    H = hamiltonian(static, [], basis=basis, dtype=np.float64, **nc).toarray()
    fu = hamiltonian([["n|", [[1.0, i] for i in range(n)]]], [], basis=basis,
                     dtype=np.float64, **nc).diagonal().real
    fd = hamiltonian([["|n", [[1.0, i] for i in range(n)]]], [], basis=basis,
                     dtype=np.float64, **nc).diagonal().real

    def cdag_up(a):
        return hamiltonian([["+|", [[1.0, a]]]], [], basis=basis, dtype=np.float64, **nc).toarray()

    def cdag_dn(a):
        return hamiltonian([["|+", [[1.0, a]]]], [], basis=basis, dtype=np.float64, **nc).toarray()

    def num_up(a):
        return hamiltonian([["n|", [[1.0, a]]]], [], basis=basis, dtype=np.float64, **nc).toarray()

    def num_dn(a):
        return hamiltonian([["|n", [[1.0, a]]]], [], basis=basis, dtype=np.float64, **nc).toarray()

    def hop_up(i, j):   # c^+_{i up} c_{j up}
        return hamiltonian([["+-|", [[1.0, i, j]]]], [], basis=basis, dtype=np.float64, **nc).toarray()

    def hop_dn(i, j):
        return hamiltonian([["|+-", [[1.0, i, j]]]], [], basis=basis, dtype=np.float64, **nc).toarray()

    def pair_op(F):     # sum_ij F_ij c^+_{i up} c^+_{j dn}
        terms = [[F[i, j], i, j] for i in range(n) for j in range(n) if abs(F[i, j]) > 1e-15]
        return hamiltonian([["+|+", terms]], [], basis=basis, dtype=np.float64, **nc).toarray()

    def cdcd(i, j):     # c^+_{i up} c^+_{j dn}
        return hamiltonian([["+|+", [[1.0, i, j]]]], [], basis=basis, dtype=np.float64, **nc).toarray()

    def cc(i, j):       # c_{i up} c_{j dn}
        return hamiltonian([["-|-", [[1.0, i, j]]]], [], basis=basis, dtype=np.float64, **nc).toarray()

    return dict(basis=basis, H=H, fu=fu, fd=fd, cdag_up=cdag_up, cdag_dn=cdag_dn,
                hop_up=hop_up, hop_dn=hop_dn, pair_op=pair_op, cdcd=cdcd, cc=cc, n=n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=2); ap.add_argument("--ly", type=int, default=2)
    ap.add_argument("--nup", type=int, default=1); ap.add_argument("--ndn", type=int, default=1)
    ap.add_argument("--U", type=float, default=4.0); ap.add_argument("--t0", type=float, default=1.0)
    ap.add_argument("--tam", type=float, default=0.0); ap.add_argument("--t1", type=float, default=0.0)
    ap.add_argument("--lam", type=float, default=0.15); ap.add_argument("--mu", type=float, default=0.0)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    assert a.nup == a.ndn, "AGP singlet trial needs nup == ndn (N pairs)"
    N = a.nup
    n = a.lx * a.ly

    Ku, Kd = am_hopping(a.lx, a.ly, a.t0, a.tam, a.t1)
    F, info = dwave_geminal(a.lx, a.ly, a.nup, a.ndn, t0=a.t0, tam=a.tam, t1=a.t1,
                            lam=a.lam, mu=a.mu)
    print(f"# {a.lx}x{a.ly} N={N}+{N} U={a.U} t1={a.t1} lam={a.lam} mu={a.mu}")
    print(f"# BdG: herm_err={info['bdg_herm_err']:.1e} min|eval|={info['min_abs_eval']:.3e} "
          f"cond(U)={info['Ucond']:.2e}")

    # a Slater walker: free-electron orbitals of Ku/Kd plus a small random rotation
    rng = np.random.default_rng(a.seed)
    wu, vu = np.linalg.eigh(Ku); wd, vd = np.linalg.eigh(Kd)
    A = vu[:, :N].copy(); B = vd[:, :N].copy()
    A = A + 0.05 * rng.standard_normal(A.shape)      # generic (non-eigen) walker
    B = B + 0.05 * rng.standard_normal(B.shape)
    A, _ = np.linalg.qr(A); B, _ = np.linalg.qr(B)

    # ---- formula ----
    c = agp_contractions(F, A, B)
    O_f, Gu_f, Gd_f, Kd_f, Ka_f = c["O"], c["Gu"], c["Gd"], c["Kd"], c["Ka"]

    # ---- exact (QuSpin) ----
    ops = build_ops(a.lx, a.ly, Ku, Kd, a.U)
    vac = np.zeros(ops["basis"].Ns)
    iv = np.where((np.rint(ops["fu"]) == 0) & (np.rint(ops["fd"]) == 0))[0]
    assert len(iv) == 1
    vac[iv[0]] = 1.0
    # |Psi_AGP> = (pair_op)^N |0>
    Pop = ops["pair_op"](F)
    psi = vac.copy()
    for _ in range(N):
        psi = Pop @ psi
    # |Phi> = prod_m c^+(A col m, up) prod_m c^+(B col m, dn) |0>
    phi = vac.copy()
    for m in range(N):
        cu = sum(A[a_, m] * ops["cdag_up"](a_) for a_ in range(n))
        phi = cu @ phi
    for m in range(N):
        cd = sum(B[b_, m] * ops["cdag_dn"](b_) for b_ in range(n))
        phi = cd @ phi

    O_ed = float(psi @ phi)
    if abs(O_ed) < 1e-300:
        print("ERROR: exact overlap ~ 0 (walker orthogonal to AGP); pick another seed.")
        return
    Gu_ed = np.array([[float(psi @ (ops["hop_up"](i, j) @ phi)) for j in range(n)] for i in range(n)]) / O_ed
    Gd_ed = np.array([[float(psi @ (ops["hop_dn"](i, j) @ phi)) for j in range(n)] for i in range(n)]) / O_ed
    E_ed = float(psi @ (ops["H"] @ phi)) / O_ed
    # anomalous tensors need the UNPROJECTED BCS bra (the projected AGP gives 0):
    # E_local = <BCS|H|Phi>/<BCS|Phi> = <Psi_AGP|H|Phi>/<Psi_AGP|Phi> (H conserves N).
    from scipy.linalg import expm
    bcs = expm(Pop) @ vac
    Ob = float(bcs @ phi)
    Kd_ed = np.array([[float(bcs @ (ops["cdcd"](i, j) @ phi)) for j in range(n)] for i in range(n)]) / Ob
    Ka_ed = np.array([[float(bcs @ (ops["cc"](i, j) @ phi)) for j in range(n)] for i in range(n)]) / Ob

    # rigorous local energy from agp.py (normal + anomalous Wick)
    model = LatticeModel(Ku, a.lx, a.ly, 0.01, a.U, 0.0, 0.0, K_dn=Kd)
    E_f = agp_local_energy(model, F, A, B)

    # ---- report ----
    print(f"\noverlap:  formula={O_f:+.6e}  exact={O_ed:+.6e}  ratio={O_f/O_ed:+.6f}"
          f"  (const incl. 1/N! and sign)")
    print(f"trace check: tr(Gu_f)={np.trace(Gu_f):.6f} tr(Gd_f)={np.trace(Gd_f):.6f}  (=N={N})")
    du = np.max(np.abs(Gu_f - Gu_ed)); dd = np.max(np.abs(Gd_f - Gd_ed))
    dKd = np.max(np.abs(Kd_f - Kd_ed)); dKa = np.max(np.abs(Ka_f - Ka_ed))
    print(f"max|Gu_formula - Gu_exact|  = {du:.3e}")
    print(f"max|Gd_formula - Gd_exact|  = {dd:.3e}")
    print(f"max|Kd_formula - Kd_exact|  = {dKd:.3e}   (anomalous <c+_u c+_d>)")
    print(f"max|Ka_formula - Ka_exact|  = {dKa:.3e}   (anomalous <c_u c_d>)")
    print(f"local energy (normal+anomalous): formula={E_f:.6f}  exact={E_ed:.6f}  |diff|={abs(E_f-E_ed):.3e}")
    ok = du < 1e-9 and dd < 1e-9 and dKd < 1e-9 and dKa < 1e-9 and abs(E_f - E_ed) < 1e-8
    print(f"\n{'PASS' if ok else 'FAIL'}: AGP overlap / Green / anomalous / energy formulas vs ED")


if __name__ == "__main__":
    main()
