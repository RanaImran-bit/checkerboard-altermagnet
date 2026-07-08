#!/usr/bin/env python3
"""CP-AFQMC with the AGP (number-projected-BCS) trial -- attack the d-wave vertex bias.

The free-electron trial has the wrong nodes for the 4-point d-wave channel: on the
4x2/t1=0.3/U=4 gate the equal-time d-wave vertex is ED +23.18 vs CP-AFQMC(fixed) ~ -2
(sign-flipped). This driver swaps in the pairing-aware AGP trial (pyqmc/agp.py,
overlap det(A^T F B), validated local energy with anomalous terms) as BOTH the
constraint/importance function AND the mixed-estimator bra, and measures:

  * energy   -- gate vs ED (validates the whole AGP-trial AFQMC end-to-end)
  * d-wave pair structure factor S_d = <Delta_d^dag Delta_d>, full + connected vertex,
    using the AGP generalized Wick (normal Gu/Gd plus anomalous Kd/Ka):
        S_d^full(walker) = sum Gu (Fd Gd Fd^T)  -  (Fd:Kd)(Fd:Ka)
    vertex = full - disconnected(ensemble-averaged Gu,Gd,Kd,Ka).

Walkers stay number-conserving Slater determinants; only the trial is paired. The
trial overlap can go negative, so this is a constrained-path (q=max(R,0)) run with
the AGP overlap as the importance function. cpqmc.py is untouched.

GEMINAL CHOICE. `--eta E` uses the WELL-CONDITIONED augmented geminal
F = phi_up phi_dn^T + E*Fd: eta=0 reproduces the free-electron trial EXACTLY (gate),
eta>0 turns on d-wave pairing smoothly. (Without --eta the bare nodal BdG geminal is
used; it has v/u->inf node directions, needs Tikhonov reg, and is a POOR importance
function -- see the report. Prefer --eta.)

STATUS (2026-06-25). Validated: U=0 -> energy exact and d-wave VERTEX = 0 exactly;
eta=0 -> reproduces the free-electron-trial CP-AFQMC energy (= ED on clean clusters).
So the machinery (AGP overlap importance + constraint + normal/anomalous estimators)
is correct. OPEN: a static eta>0 importance function is variance-unstable at half
filling and does not yet cleanly reduce the d-wave-vertex CP bias; next steps are a
back-propagated (not mixed) AGP estimator, eta optimization / variance control, and a
full-Fock-space ED pairing reference (cf. validate_chid) for an absolute comparison.

    source tools/env.sh
    # gate: eta=0 energy must match ED, U=0 vertex must be 0
    python pyqmc/agp_dwave.py --lx 2 --ly 2 --nup 1 --ndn 1 --U 4 --t1 0.3 --eta 0.0 --ed
    # turn on pairing
    python pyqmc/agp_dwave.py --lx 4 --ly 2 --nup 4 --ndn 4 --U 4 --t1 0.3 --eta 0.2
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
import argparse, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from cpqmc import LatticeModel, am_hopping, _green
from agp import dwave_geminal, augmented_geminal, agp_contractions, agp_local_energy


def build_Fd(lx, ly):
    """d-wave bond form-factor matrix Fd[m,j] (+x bonds, -y bonds); s-wave Fs."""
    n = lx * ly
    Fs = np.zeros((n, n)); Fd = np.zeros((n, n))
    def idx(x, y): return (x % lx) * ly + (y % ly)
    for x in range(lx):
        for y in range(ly):
            m = idx(x, y)
            for (dx, dy, fd) in ((1, 0, 1.0), (-1, 0, 1.0), (0, 1, -1.0), (0, -1, -1.0)):
                j = idx(x + dx, y + dy); Fs[m, j] += 1.0; Fd[m, j] += fd
    return Fs, Fd


class AGPTrial:
    """Thin AGP-trial wrapper exposing the overlap (for the CP constraint/importance)
    and the contraction set (for estimators)."""
    def __init__(self, F):
        self.F = F
    def overlap(self, A, B):
        M = A.T @ self.F @ B
        s, ld = np.linalg.slogdet(M)
        return s * np.exp(ld)


def step_pop(model, trial, pu, pd, w, rng):
    """One constrained-path propagation step over the population with the AGP
    overlap as the importance function (q=max(R,0), positive weights)."""
    eu, ed = model.expK, model.expK_dn
    for i in range(len(w)):
        if w[i] == 0:
            continue
        u = eu @ pu[i]; d = ed @ pd[i]
        O = trial.overlap(u, d)
        if O == 0:
            w[i] = 0.0; continue
        dead = False
        for (a, sa, b, sb, e1, e2, cf) in model.terms:
            R = np.empty(2); cand = []
            for k in (0, 1):
                u2 = u.copy(); d2 = d.copy()
                (u2 if sa == 0 else d2)[a, :] *= e1[k]
                (u2 if sb == 0 else d2)[b, :] *= e2[k]
                O2 = trial.overlap(u2, d2)
                R[k] = cf[k] * O2 / O
                cand.append((u2, d2, O2))
            q = np.maximum(R, 0.0)
            tot = q.sum()
            if tot <= 0:
                w[i] = 0.0; dead = True; break
            k = rng.choice(2, p=q / tot)
            u, d, O = cand[k]
            w[i] *= 0.5 * tot
        if dead:
            continue
        pu[i] = eu @ u; pd[i] = ed @ d
    return pu, pd, w


def measure(model, trial, pu, pd, w, Fd):
    """Mixed-estimator energy and the equal-time d-wave pair structure factor
    S_d = <Delta_d Delta_d^dag> (full + connected vertex) with the AGP bra, in the
    SAME ordering as the ED reference / validate_chid C(tau=0). Generalized Wick:
        S_full = sum Pu (Fd Pd Fd^T)  -  (Fd:Ka)(Fd:Kd) ,   Pu = I - Gu^T  (particle GF)
    vertex = S_full - bubble, bubble = sum Pu_avg (Fd Pd_avg Fd^T) (normal disconnected;
    the anomalous -(Fd:Ka)(Fd:Kd) is a connected pairing contribution and stays in)."""
    n = model.n
    F = trial.F
    I = np.eye(n)
    sw = 0.0; E = 0.0; Sfull = 0.0
    Pu_a = np.zeros((n, n)); Pd_a = np.zeros((n, n))
    amp_dag = 0.0; amp = 0.0
    for i in range(len(w)):
        if w[i] == 0:
            continue
        try:
            c = agp_contractions(F, pu[i], pd[i])
        except np.linalg.LinAlgError:
            continue
        Gu, Gd, Kd, Ka = c["Gu"], c["Gd"], c["Kd"], c["Ka"]
        wi = w[i]
        E += wi * agp_local_energy(model, F, pu[i], pd[i])
        Pu = I - Gu.T; Pd = I - Gd.T            # particle GFs <c c^+>
        normal = float((Pu * (Fd @ Pd @ Fd.T)).sum())
        FKd = float((Fd * Kd).sum()); FKa = float((Fd * Ka).sum())
        Sfull += wi * (normal - FKa * FKd)
        Pu_a += wi * Pu; Pd_a += wi * Pd
        amp_dag += wi * FKd; amp += wi * FKa
        sw += wi
    if sw == 0:
        return dict(nan=True)
    E /= sw; Sfull /= sw; Pu_a /= sw; Pd_a /= sw; amp_dag /= sw; amp /= sw
    bubble = float((Pu_a * (Fd @ Pd_a @ Fd.T)).sum())   # normal disconnected
    return dict(E=E, Sd_full=Sfull, Sd_vtx=Sfull - bubble,
                pair_amp=amp_dag, nan=False)


def one_pop(P, seed):
    Ku, Kd = am_hopping(P["lx"], P["ly"], P["t0"], P["tam"], P["t1"])
    model = LatticeModel(Ku, P["lx"], P["ly"], P["dt"], P["U"], 0.0, 0.0, K_dn=Kd)
    _, Fd = build_Fd(P["lx"], P["ly"])
    vu, vd = model.eigvecs, model.eigvecs_dn
    A0 = vu[:, :P["nup"]].copy(); B0 = vd[:, :P["ndn"]].copy()
    if P["eta"] is None:                          # bare BdG geminal (nodal; needs reg)
        F, info = dwave_geminal(P["lx"], P["ly"], P["nup"], P["ndn"], t0=P["t0"],
                                tam=P["tam"], t1=P["t1"], lam=P["lam"], mu=P["mu"])
        Fnorm = info["Fnorm"]
    else:                                         # augmented: free trial + eta*Fd
        F = augmented_geminal(A0, B0, Fd, P["eta"]); Fnorm = float(np.linalg.norm(F))
    trial = AGPTrial(F)
    rng = np.random.default_rng(seed)
    nw = P["nw"]
    # start walkers from the free-electron determinant (a valid fixed-N ket)
    pu = np.stack([A0.copy() for _ in range(nw)])
    pd = np.stack([B0.copy() for _ in range(nw)])
    w = np.ones(nw)
    for it in range(P["nequil"] + P["nmeas"]):
        pu, pd, w = step_pop(model, trial, pu, pd, w, rng)
        if (it + 1) % 10 == 0:                      # reorth + comb pop control
            for i in range(nw):
                if w[i] > 0:
                    pu[i] = np.linalg.qr(pu[i])[0]; pd[i] = np.linalg.qr(pd[i])[0]
            wc = np.clip(w, 0, None); s = wc.sum()
            if s > 0:
                wc *= nw / s; cum = np.cumsum(wc) / wc.sum()
                idx = np.clip(np.searchsorted(cum, (rng.random() / nw) + np.arange(nw) / nw), 0, nw - 1)
                pu = pu[idx].copy(); pd = pd[idx].copy(); w = np.ones(nw)
    out = measure(model, trial, pu, pd, w, Fd)
    out["Fnorm"] = Fnorm
    return out


def ed_reference(P):
    """FULL-FOCK-SPACE ED gate (SPARSE): ground-state energy + the EQUAL-TIME d-wave
    pair structure factor S_d = <Delta_d Delta_d^dag> (full) and its connected VERTEX,
    with the SPIN-DEPENDENT am_hopping (Ku/Kd, incl. t1). Delta_d^dag raises N so this
    needs the full Fock space; sparse so it scales past 2x2 (dense 4^8 = 34 GB).
        S_full = || Delta_d^dag |0> ||^2 ,  Gu0_ab = <c_a up c^+_b up> (particle GF),
        bubble = sum Gu0 (Fd Gd0 Fd^T) ,  vertex = S_full - bubble.
    The (nup,ndn) ground state is selected by a number penalty
    lam*((N_up-nup)^2 + (N_dn-ndn)^2) added to H; Gu0 is the Gram matrix of the n
    vectors c^+_{b up}|0> (so only 2n operator applications, no dense n^2 ops).
    Returns (E0, S_full, S_vertex)."""
    import scipy.sparse as sp
    from scipy.sparse.linalg import eigsh
    from quspin.basis import spinful_fermion_basis_general
    from quspin.operators import hamiltonian
    lx, ly = P["lx"], P["ly"]; n = lx * ly
    nup, ndn, U = P["nup"], P["ndn"], P["U"]
    Ku, Kd = am_hopping(lx, ly, P["t0"], P["tam"], P["t1"])
    basis = spinful_fermion_basis_general(n)              # full Fock space
    nc = dict(check_pcon=False, check_symm=False, check_herm=False)
    up = [[Ku[i, j], i, j] for i in range(n) for j in range(n) if abs(Ku[i, j]) > 1e-15]
    dn = [[Kd[i, j], i, j] for i in range(n) for j in range(n) if abs(Kd[i, j]) > 1e-15]
    inter = [[U, i, i] for i in range(n)]
    H = hamiltonian([["+-|", up], ["|+-", dn], ["n|n", inter]], [], basis=basis,
                    dtype=np.float64, **nc).tocsr()
    fu = hamiltonian([["n|", [[1.0, i] for i in range(n)]]], [], basis=basis,
                     dtype=np.float64, **nc).diagonal().real
    fd = hamiltonian([["|n", [[1.0, i] for i in range(n)]]], [], basis=basis,
                     dtype=np.float64, **nc).diagonal().real
    lam = 100.0
    pen = sp.diags(lam * ((fu - nup) ** 2 + (fd - ndn) ** 2))
    E, V = eigsh((H + pen).tocsc(), k=1, which="SA")
    psi = V[:, 0]; E0 = float(psi @ (H @ psi))
    _, Fd = build_Fd(lx, ly)
    bonds = [(m, j, Fd[m, j]) for m in range(n) for j in range(n) if abs(Fd[m, j]) > 1e-15]
    Dd = hamiltonian([["+|+", [[f, m, j] for (m, j, f) in bonds]]], [], basis=basis,
                     dtype=np.float64, **nc).tocsr()
    dpsi = Dd @ psi                                       # Delta_d^dag |0>  (raises N)
    S_full = float(dpsi @ dpsi)
    # particle GFs at tau=0 via the Gram trick: Gu0_ab = <c^+_a up psi | c^+_b up psi>
    Vu = np.empty((basis.Ns, n)); Vd = np.empty((basis.Ns, n))
    for b in range(n):
        Vu[:, b] = hamiltonian([["+|", [[1.0, b]]]], [], basis=basis, dtype=np.float64, **nc).tocsr() @ psi
        Vd[:, b] = hamiltonian([["|+", [[1.0, b]]]], [], basis=basis, dtype=np.float64, **nc).tocsr() @ psi
    Gu0 = Vu.T @ Vu; Gd0 = Vd.T @ Vd
    bubble = float((Gu0 * (Fd @ Gd0 @ Fd.T)).sum())
    return E0, S_full, S_full - bubble


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=4); ap.add_argument("--ly", type=int, default=2)
    ap.add_argument("--nup", type=int, default=4); ap.add_argument("--ndn", type=int, default=4)
    ap.add_argument("--U", type=float, default=4.0); ap.add_argument("--t0", type=float, default=1.0)
    ap.add_argument("--tam", type=float, default=0.0); ap.add_argument("--t1", type=float, default=0.3)
    ap.add_argument("--lam", type=float, default=0.6); ap.add_argument("--mu", type=float, default=0.0)
    ap.add_argument("--eta", type=float, default=None,
                    help="augmented geminal F=phi_up phi_dn^T + eta*Fd (eta=0 -> free trial). "
                         "If unset, use the bare (nodal, regularized) BdG geminal.")
    ap.add_argument("--dt", type=float, default=0.02); ap.add_argument("--nw", type=int, default=80)
    ap.add_argument("--nequil", type=int, default=120); ap.add_argument("--nmeas", type=int, default=1)
    ap.add_argument("--npop", type=int, default=24); ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--ed", action="store_true")
    a = ap.parse_args()
    assert a.nup == a.ndn
    P = dict(lx=a.lx, ly=a.ly, nup=a.nup, ndn=a.ndn, U=a.U, t0=a.t0, tam=a.tam, t1=a.t1,
             lam=a.lam, mu=a.mu, eta=a.eta, dt=a.dt, nw=a.nw, nequil=a.nequil, nmeas=a.nmeas)
    import multiprocessing as mp
    with mp.Pool(min(a.npop, mp.cpu_count())) as pool:
        res = pool.starmap(one_pop, [(P, a.seed + p) for p in range(a.npop)])
    res = [r for r in res if not r.get("nan")]
    def stat(key):
        v = np.array([r[key] for r in res]); return v.mean(), v.std() / np.sqrt(len(v))
    print(f"# {a.lx}x{a.ly} nup={a.nup} ndn={a.ndn} U={a.U} tam={a.tam} t1={a.t1} "
          f"lam={a.lam} mu={a.mu} dt={a.dt} npop={a.npop} nw={a.nw}")
    print(f"# AGP-trial CP-AFQMC. |F|={res[0]['Fnorm']:.3g}")
    for key, name in (("E", "energy"), ("Sd_full", "S_d full"), ("Sd_vtx", "S_d VERTEX"),
                      ("pair_amp", "<Delta_d^dag>")):
        m, e = stat(key); print(f"{name:14s} = {m:+.5f} +/- {e:.5f}")
    if a.ed:
        E0, Sf, Sv = ed_reference(P)
        print(f"\n# full-Fock ED gate:")
        print(f"{'energy':14s} = {E0:+.5f}   (eta=0 AGP energy must match)")
        print(f"{'S_d full':14s} = {Sf:+.5f}")
        print(f"{'S_d VERTEX':14s} = {Sv:+.5f}")


if __name__ == "__main__":
    main()
