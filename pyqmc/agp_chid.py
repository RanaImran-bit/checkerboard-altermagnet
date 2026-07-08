#!/usr/bin/env python3
"""BACK-PROPAGATED AGP-bra d-wave PAIRING SUSCEPTIBILITY chi_d = integral C_d(tau) dtau.

The AGP-bra analogue of cpqmc.run_bp_chid: the unequal-time singlet pairing correlator
C_d(tau) = <Delta_d(tau) Delta_d^dag(0)> measured with a back-propagated AGP (number-
projected-BCS) bra instead of the free-electron bra. Generalized Wick (paired bra,
Slater ket) gives, per walker,

    C_d(tau_l) = sum Pu^tau (Fd Pd^tau Fd^T)  -  (Fd:Ka^tau)(Fd:Kd^0)
                  \_______ normal kernel _______/   \__ anomalous (paired bra) __/

with the time-displaced particle GFs  Pu^tau = B_(l)(I - Gu^T),  Pd^tau = B_(l)(I - Gd^T)
and the time-displaced equal-time anomalous GF  Ka^tau = Bu_(l) Ka^0 Bd_(l)^T.  Gu,Gd,
Kd,Ka are the AGP contractions at slice 0 from the FULL-window back-propagated geminal
F_bp = Bu_tot^T F Bd_tot (exactly as agp_bp_vertex). At tau=0 this reduces to the
ED-validated equal-time d-wave vertex. The connected VERTEX subtracts the disconnected
normal bubble AND the disconnected anomalous product (built from ensemble-averaged GFs).

CONSTRAINT/importance = stable free-electron determinant (AGP-as-constraint is unstable);
the AGP enters only as the measurement bra. eta=0 reproduces the free-bra chid (gate).

    source tools/env.sh
    python pyqmc/agp_chid.py --lx 2 --ly 2 --nup 1 --ndn 1 --U 0 --etas 0.0,0.5   # U=0 gate
    python pyqmc/agp_chid.py --lx 2 --ly 2 --nup 1 --ndn 1 --U 4 --etas 0.0,0.5 --ed
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
import argparse, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from cpqmc import LatticeModel, am_hopping
from agp import augmented_geminal, agp_contractions
from agp_dwave import build_Fd
from agp_bp_vertex import step_record_free, comb


def step_bmats(model, ch):
    """Per-spin one-body forward propagators bu=e^{-dtK/2}D_up e^{-dtK/2} (dn likewise)
    for one recorded HS-field index list ch (model.terms order). Mirrors
    cpqmc.Estimators._step_bmats."""
    n = model.n
    d_up = np.ones(n); d_dn = np.ones(n)
    for (a, sa, b, sb, e1, e2, cf), k in zip(model.terms, ch):
        (d_up if sa == 0 else d_dn)[a] *= e1[k]
        (d_up if sb == 0 else d_dn)[b] *= e2[k]
    bu = model.expK @ (d_up[:, None] * model.expK)
    bd = model.expK_dn @ (d_dn[:, None] * model.expK_dn)
    return bu, bd


def chid_block_agp(model, Phi_u, Phi_d, ket_up, ket_dn, rec, w, Fd, etas, bp):
    """Per-block weighted C_d(tau_l), l=0..bp, for each eta (F = Phi_u Phi_d^T + eta Fd).
    Returns dict eta -> (Cfull[L], Cbub[L], W) accumulators (full and disconnected
    kernels). Disconnected built from ensemble-averaged Pu^tau, Pd^tau, FKa^tau, FKd^0."""
    n = model.n; L = bp + 1; I = np.eye(n)
    out = {}
    for eta in etas:
        F = augmented_geminal(Phi_u, Phi_d, Fd, eta)
        Cfull = np.zeros(L); W = 0.0
        Pu_avg = np.zeros((L, n, n)); Pd_avg = np.zeros((L, n, n))
        FKa_avg = np.zeros(L); FKd_avg = 0.0
        for i in range(len(w)):
            if w[i] == 0 or any(c is None for c in rec[i]):
                continue
            # full-window back-propagated geminal -> AGP contractions at slice 0
            Bu_tot = np.eye(n); Bd_tot = np.eye(n)
            bmats = []
            for ch in rec[i]:
                bu, bd = step_bmats(model, ch); bmats.append((bu, bd))
                Bu_tot = bu @ Bu_tot; Bd_tot = bd @ Bd_tot
            F_bp = Bu_tot.T @ F @ Bd_tot
            try:
                c = agp_contractions(F_bp, ket_up[i], ket_dn[i])
            except np.linalg.LinAlgError:
                continue
            Gu, Gd, Kd, Ka = c["Gu"], c["Gd"], c["Kd"], c["Ka"]
            Pu0 = I - Gu.T; Pd0 = I - Gd.T
            FKd0 = float((Fd * Kd).sum())
            wi = w[i]
            Bu = np.eye(n); Bd = np.eye(n)
            for l in range(L):
                if l > 0:
                    bu, bd = bmats[l - 1]
                    Bu = bu @ Bu; Bd = bd @ Bd
                Pu = Bu @ Pu0; Pd = Bd @ Pd0
                Kat = Bu @ Ka @ Bd.T
                FKat = float((Fd * Kat).sum())
                Cfull[l] += wi * (float((Pu * (Fd @ Pd @ Fd.T)).sum()) - FKat * FKd0)
                Pu_avg[l] += wi * Pu; Pd_avg[l] += wi * Pd; FKa_avg[l] += wi * FKat
            FKd_avg += wi * FKd0
            W += wi
        out[eta] = (Cfull, Pu_avg, Pd_avg, FKa_avg, FKd_avg, W)
    return out


def one_pop(P, seed):
    Ku, Kd = am_hopping(P["lx"], P["ly"], P["t0"], P["tam"], P["t1"])
    model = LatticeModel(Ku, P["lx"], P["ly"], P["dt"], P["U"], 0.0, 0.0, K_dn=Kd)
    Phi_u = model.eigvecs[:, :P["nup"]].copy(); Phi_d = model.eigvecs_dn[:, :P["ndn"]].copy()
    _, Fd = build_Fd(P["lx"], P["ly"])
    etas = P["etas"]; bp = P["bp"]; L = bp + 1
    rng = np.random.default_rng(seed)
    nw = P["nw"]
    pu = np.stack([Phi_u.copy() for _ in range(nw)]); pd = np.stack([Phi_d.copy() for _ in range(nw)])
    w = np.ones(nw)
    for _ in range(P["nequil"]):
        rec0 = [[] for _ in range(nw)]
        pu, pd, w = step_record_free(model, Phi_u, Phi_d, pu, pd, w, rng, rec0)
        pu, pd, w = comb(model, pu, pd, w, rng)
    tot = {e: [np.zeros(L), np.zeros((L, model.n, model.n)), np.zeros((L, model.n, model.n)),
               np.zeros(L), 0.0, 0.0] for e in etas}
    for blk in range(P["nblocks"]):
        ket_up = pu.copy(); ket_dn = pd.copy()
        rec = [[] for _ in range(nw)]
        for _ in range(bp):
            pu, pd, w = step_record_free(model, Phi_u, Phi_d, pu, pd, w, rng, rec)
        res = chid_block_agp(model, Phi_u, Phi_d, ket_up, ket_dn, rec, w, Fd, etas, bp)
        for e in etas:
            Cf, Pua, Pda, FKaa, FKd, W = res[e]
            if W > 0:
                tot[e][0] += Cf; tot[e][1] += Pua; tot[e][2] += Pda
                tot[e][3] += FKaa; tot[e][4] += FKd; tot[e][5] += W
        pu, pd, w = comb(model, pu, pd, w, rng)
    # per-eta: full + vertex C(tau), then chi = trapezoid
    dt = P["dt"]; taus = dt * np.arange(L)
    res = {}
    for e in etas:
        Cf, Pua, Pda, FKaa, FKd, W = tot[e]
        if W <= 0:
            res[e] = (np.full(L, np.nan), np.full(L, np.nan)); continue
        full = Cf / W
        Pua /= W; Pda /= W; FKaa /= W; FKd /= W
        bub = np.array([float((Pua[l] * (Fd @ Pda[l] @ Fd.T)).sum()) - FKaa[l] * FKd
                        for l in range(L)])
        res[e] = (full, full - bub)            # full, vertex C(tau)
    return res, taus


def trapz(C, dt):
    return dt * (C[1:-1].sum() + 0.5 * (C[0] + C[-1]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=2); ap.add_argument("--ly", type=int, default=2)
    ap.add_argument("--nup", type=int, default=1); ap.add_argument("--ndn", type=int, default=1)
    ap.add_argument("--U", type=float, default=4.0); ap.add_argument("--t0", type=float, default=1.0)
    ap.add_argument("--tam", type=float, default=0.0); ap.add_argument("--t1", type=float, default=0.0)
    ap.add_argument("--etas", default="0.0,0.5")
    ap.add_argument("--dt", type=float, default=0.05); ap.add_argument("--nw", type=int, default=120)
    ap.add_argument("--nequil", type=int, default=120); ap.add_argument("--nblocks", type=int, default=40)
    ap.add_argument("--bp", type=int, default=16); ap.add_argument("--npop", type=int, default=12)
    ap.add_argument("--seed", type=int, default=1); ap.add_argument("--ed", action="store_true")
    a = ap.parse_args()
    etas = [float(x) for x in a.etas.split(",")]
    P = dict(lx=a.lx, ly=a.ly, nup=a.nup, ndn=a.ndn, U=a.U, t0=a.t0, tam=a.tam, t1=a.t1,
             etas=etas, dt=a.dt, nw=a.nw, nequil=a.nequil, nblocks=a.nblocks, bp=a.bp)
    import multiprocessing as mp
    with mp.Pool(min(a.npop, mp.cpu_count())) as pool:
        outs = pool.starmap(one_pop, [(P, a.seed + p) for p in range(a.npop)])
    taus = outs[0][1]
    print(f"# {a.lx}x{a.ly} nup={a.nup} ndn={a.ndn} U={a.U} tam={a.tam} t1={a.t1} dt={a.dt} "
          f"bp={a.bp} npop={a.npop}")
    ed = None
    if a.ed:
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "ed"))
        from validate_chid import ed_chid
        _, edres = ed_chid(a.lx, a.ly, a.nup, a.ndn, a.t0, a.U, taus)
        ed = (trapz(edres["d"][0], a.dt), trapz(edres["d"][2], a.dt))  # full, vertex
        print(f"# ED (full-Fock Lehmann): chi_d FULL={ed[0]:.4f}  chi_d VERTEX={ed[1]:.4f}")
    print(f"{'eta':>6} {'chi_d FULL':>12} {'+/-':>8} {'chi_d VERTEX':>13} {'+/-':>8}")
    for e in etas:
        fulls = np.array([trapz(o[0][e][0], a.dt) for o in outs])
        vtxs = np.array([trapz(o[0][e][1], a.dt) for o in outs])
        def stat(x):
            x = x[np.isfinite(x)]; return x.mean(), (x.std() / np.sqrt(len(x)) if len(x) > 1 else np.nan)
        fm, fe = stat(fulls); vm, ve = stat(vtxs)
        print(f"{e:>6.2f} {fm:>12.4f} {fe:>8.4f} {vm:>13.4f} {ve:>8.4f}")


if __name__ == "__main__":
    main()
