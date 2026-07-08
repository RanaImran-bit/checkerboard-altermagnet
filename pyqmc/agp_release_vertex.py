#!/usr/bin/env python3
"""Constraint RELEASE (ket) + fixed AGP/BCS bra -> close the d-wave vertex KET bias.

Pure back-propagation with the AGP bra (agp_bp_vertex.py) saturates at ~80% of the ED
d-wave vertex: the bra is fixed but the KET is still free-trial-constrained (the CP
ceiling). This driver removes that last piece by RELEASING the constraint on the ket,
exactly as pyqmc/release_dwave.py (free signed projection after CP equilibration), but
measuring the equal-time d-wave vertex with the pairing-aware AGP bra instead of the
free-electron bra:

    vertex(tau) = <Psi_AGP| Delta_d Delta_d^dag |phi(tau)> / <Psi_AGP|phi(tau)>  (connected)

As release time tau grows the released walker phi(tau) -> exact GS (ket bias removed),
while the AGP bra carries the pairing channel (anomalous Kd/Ka) the free bra lacks --
so the vertex should climb past the CP ceiling toward ED, until <sign> decays. tau=0 is
the CP value; eta=0 is the free-bra control. cpqmc.py untouched.

    source tools/env.sh
    python pyqmc/agp_release_vertex.py --lx 2 --ly 2 --nup 2 --ndn 2 --U 4 --t1 0.3 --eta 0.5 --ed
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
import argparse, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from cpqmc import LatticeModel, TrialWF, am_hopping
from agp import augmented_geminal, agp_contractions
from agp_dwave import build_Fd, ed_reference


def step_pop(model, trial, pu, pd, w, rng, constrained):
    """One propagation step over the population using the FREE-electron trial overlap
    as the importance function. constrained: q=max(R,0) (CP); else q=|R| signed
    (free projection / release)."""
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
            q = np.maximum(R, 0.0) if constrained else np.abs(R)
            tot = q.sum()
            if tot <= 0:
                w[i] = 0.0; dead = True; break
            k = rng.choice(2, p=q / tot)
            u, d, O = cand[k]
            w[i] *= (0.5 * tot) if constrained else (0.5 * np.sign(R[k]) * np.abs(R).sum())
        if dead:
            continue
        pu[i] = eu @ u; pd[i] = ed @ d
    return pu, pd, w


def measure_agp_vertex(F, pu, pd, w, Fd):
    """Connected equal-time d-wave vertex with the fixed AGP bra against the current
    (released) walkers. Signed-weighted; returns (full, vertex, sign)."""
    n = Fd.shape[0]; I = np.eye(n)
    sw = 0.0; saw = 0.0; Sfull = 0.0
    Pu_a = np.zeros((n, n)); Pd_a = np.zeros((n, n)); FKd_a = 0.0; FKa_a = 0.0
    for i in range(len(w)):
        if w[i] == 0:
            continue
        try:
            c = agp_contractions(F, pu[i], pd[i])
        except np.linalg.LinAlgError:
            continue
        Gu, Gd, Kd, Ka = c["Gu"], c["Gd"], c["Kd"], c["Ka"]
        Pu = I - Gu.T; Pd = I - Gd.T
        normal = float((Pu * (Fd @ Pd @ Fd.T)).sum())
        FKd = float((Fd * Kd).sum()); FKa = float((Fd * Ka).sum())
        Sfull += w[i] * (normal - FKa * FKd)
        Pu_a += w[i] * Pu; Pd_a += w[i] * Pd; FKd_a += w[i] * FKd; FKa_a += w[i] * FKa
        sw += w[i]; saw += abs(w[i])
    if abs(sw) < 1e-300:
        return np.nan, np.nan, 0.0
    Sfull /= sw; Pu_a /= sw; Pd_a /= sw; FKd_a /= sw; FKa_a /= sw
    bubble = float((Pu_a * (Fd @ Pd_a @ Fd.T)).sum())
    vtx = Sfull - bubble + FKd_a * FKa_a       # subtract normal + anomalous disconnected
    return Sfull, vtx, sw / saw


def one_pop(P, seed):
    Ku, Kd = am_hopping(P["lx"], P["ly"], P["t0"], P["tam"], P["t1"])
    model = LatticeModel(Ku, P["lx"], P["ly"], P["dt"], P["U"], 0.0, 0.0, K_dn=Kd)
    vu, vd = model.eigvecs, model.eigvecs_dn
    Phi_u = vu[:, :P["nup"]].copy(); Phi_d = vd[:, :P["ndn"]].copy()
    trial = TrialWF(Phi_u.copy(), Phi_d.copy(), mode="fixed")   # free-electron CONSTRAINT
    _, Fd = build_Fd(P["lx"], P["ly"])
    F = augmented_geminal(Phi_u, Phi_d, Fd, P["eta"])           # fixed AGP/BCS BRA
    rng = np.random.default_rng(seed)
    nw = P["nw"]
    pu = np.stack([Phi_u.copy() for _ in range(nw)])
    pd = np.stack([Phi_d.copy() for _ in range(nw)])
    w = np.ones(nw)
    for it in range(P["nequil"]):                              # CP equilibration
        pu, pd, w = step_pop(model, trial, pu, pd, w, rng, constrained=True)
        if (it + 1) % 10 == 0:
            for i in range(nw):
                if w[i] > 0:
                    pu[i] = np.linalg.qr(pu[i])[0]; pd[i] = np.linalg.qr(pd[i])[0]
            wc = np.clip(w, 0, None); s = wc.sum()
            if s > 0:
                wc *= nw / s; cum = np.cumsum(wc) / wc.sum()
                idx = np.clip(np.searchsorted(cum, (rng.random() / nw) + np.arange(nw) / nw), 0, nw - 1)
                pu = pu[idx].copy(); pd = pd[idx].copy(); w = np.ones(nw)
    out = [measure_agp_vertex(F, pu, pd, w, Fd)]               # tau=0 = CP value
    for it in range(P["nrelease"]):                            # RELEASE the ket
        pu, pd, w = step_pop(model, trial, pu, pd, w, rng, constrained=False)
        mag = np.mean(np.abs(w[w != 0])) if np.any(w != 0) else 1.0
        if mag > 0:
            w = w / mag                                        # rescale (cancels in ratios)
        out.append(measure_agp_vertex(F, pu, pd, w, Fd))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=2); ap.add_argument("--ly", type=int, default=2)
    ap.add_argument("--nup", type=int, default=2); ap.add_argument("--ndn", type=int, default=2)
    ap.add_argument("--U", type=float, default=4.0); ap.add_argument("--t0", type=float, default=1.0)
    ap.add_argument("--tam", type=float, default=0.0); ap.add_argument("--t1", type=float, default=0.3)
    ap.add_argument("--eta", type=float, default=0.5)
    ap.add_argument("--dt", type=float, default=0.02); ap.add_argument("--nw", type=int, default=80)
    ap.add_argument("--nequil", type=int, default=150); ap.add_argument("--nrelease", type=int, default=24)
    ap.add_argument("--npop", type=int, default=24); ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--ed", action="store_true")
    a = ap.parse_args()
    assert a.nup == a.ndn
    P = dict(lx=a.lx, ly=a.ly, nup=a.nup, ndn=a.ndn, U=a.U, t0=a.t0, tam=a.tam, t1=a.t1,
             eta=a.eta, dt=a.dt, nw=a.nw, nequil=a.nequil, nrelease=a.nrelease)
    import multiprocessing as mp
    with mp.Pool(min(a.npop, mp.cpu_count())) as pool:
        res = pool.starmap(one_pop, [(P, a.seed + p) for p in range(a.npop)])
    L = a.nrelease + 1
    print(f"# {a.lx}x{a.ly} nup={a.nup} ndn={a.ndn} U={a.U} tam={a.tam} t1={a.t1} eta={a.eta} "
          f"dt={a.dt} npop={a.npop} nw={a.nw}")
    print(f"# free-electron CONSTRAINT + RELEASED ket + fixed AGP/BCS bra. tau=0 = CP value.")
    if a.ed:
        E0, Sf, Sv = ed_reference(P)
        print(f"# full-Fock ED: E0={E0:.5f}  S_d full={Sf:.4f}  S_d VERTEX={Sv:.4f}")
    print(f"{'tau':>6} {'<sign>':>8} {'S_d full':>10} {'d-VERTEX':>10} {'+/-':>8}")
    for l in range(L):
        vts = np.array([r[l][1] for r in res]); fds = np.array([r[l][0] for r in res])
        sgn = np.array([r[l][2] for r in res])
        vts = vts[np.isfinite(vts)]
        m = vts.mean(); e = vts.std() / np.sqrt(len(vts)) if len(vts) > 1 else np.nan
        print(f"{l * a.dt:>6.3f} {np.nanmean(sgn):>8.4f} {np.nanmean(fds):>10.3f} {m:>10.3f} {e:>8.3f}")


if __name__ == "__main__":
    main()
