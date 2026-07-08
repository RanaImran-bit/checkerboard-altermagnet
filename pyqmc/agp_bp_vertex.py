#!/usr/bin/env python3
"""BACK-PROPAGATED AGP-bra d-wave vertex (Phase 2+3 of docs/agp_vertex_plan.md).

Decouples the constraint from the measurement bra:
  * CONSTRAINT / importance  -> the STABLE free-electron determinant (the validated,
    well-behaved CP-AFQMC ket distribution; the AGP-as-constraint is variance-unstable
    at half filling).
  * MEASUREMENT bra          -> a BACK-PROPAGATED AGP (number-projected-BCS) trial,
    which carries the anomalous pairing contractions the free bra lacks.

This isolates the "better bra" benefit: does a pairing bra alone move the equal-time
d-wave vertex toward ED, on the free-trial ket distribution, where the free bra gives
~0 / sign-flipped?

Back-propagation of the AGP bra. With the recorded forward one-body propagators
B_l = e^{-dtK/2} D_l e^{-dtK/2}, the ket at slice 0 propagates to slice L as
phi_L = B_L...B_1 phi_0. The overlap with the AGP bra obeys
  <AGP(F)| (prod B) |phi_0> = det(phi_0^T F_bp phi_0)  with  F_bp = Bu_tot^T F Bd_tot,
Bu_tot = bu_L...bu_1 (dn likewise). So the BP-AGP mixed contractions at slice 0 are
just agp_contractions(F_bp, ket_up0, ket_dn0) -- reusing the ED-validated AGP machinery.

d-wave vertex (same <Delta_d Delta_d^dag> ordering as the ED gate / validate_chid):
  S_full = sum Pu (Fd Pd Fd^T) - (Fd:Ka)(Fd:Kd),  Pu=I-Gu^T  (particle GF)
  vertex = S_full - bubble,  bubble = sum Pu_avg (Fd Pd_avg Fd^T).
Measured at several eta (F = phi_up phi_dn^T + eta Fd) from the SAME walkers/fields;
eta=0 is the free-bra control. cpqmc.py untouched.

    source tools/env.sh
    # gate: U=0 -> vertex 0 at every eta; eta=0 energy_bp matches ED
    python pyqmc/agp_bp_vertex.py --lx 2 --ly 2 --nup 2 --ndn 2 --U 4 --t1 0.3 --ed
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
import argparse, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from cpqmc import LatticeModel, am_hopping, _ov_spin
from agp import augmented_geminal, agp_contractions, agp_local_energy
from agp_dwave import build_Fd, ed_reference


def step_record_free(model, Pu_t, Pd_t, pu, pd, w, rng, rec, release=False):
    """One free-electron-importance propagation step, recording the chosen HS field
    index per term for back-propagation. Mirrors cpqmc.Propagator.step.
    release=False -> constrained path (q=max(R,0), positive weight).
    release=True  -> FREE projection (signed weight w*=0.5*sign(R_k)*sum|R|, sample
    k prop |R|) -- removes the constraint/ket bias at the cost of a decaying sign."""
    eu, ed = model.expK, model.expK_dn
    for i in range(len(w)):
        if w[i] == 0:
            rec[i].append(None); continue
        u = eu @ pu[i]; d = ed @ pd[i]
        ch = []; dead = False
        for (a, sa, b, sb, e1, e2, cf) in model.terms:
            R = np.empty(2); cand = []
            ov0 = _ov_spin(Pu_t, u) * _ov_spin(Pd_t, d)
            if ov0 == 0:
                w[i] = 0.0; dead = True; break
            for k in (0, 1):
                u2 = u.copy(); d2 = d.copy()
                (u2 if sa == 0 else d2)[a, :] *= e1[k]
                (u2 if sb == 0 else d2)[b, :] *= e2[k]
                R[k] = cf[k] * _ov_spin(Pu_t, u2) * _ov_spin(Pd_t, d2) / ov0
                cand.append((u2, d2))
            q = np.abs(R) if release else np.maximum(R, 0.0)
            tot = q.sum()
            if tot <= 0:
                w[i] = 0.0; dead = True; break
            k = int(rng.choice(2, p=q / tot))
            u, d = cand[k]
            w[i] *= (0.5 * np.sign(R[k]) * np.abs(R).sum()) if release else (0.5 * tot)
            ch.append(k)
        if dead:
            rec[i].append(None); continue
        pu[i] = eu @ u; pd[i] = ed @ d
        rec[i].append(ch)
    return pu, pd, w


def _step_bmats(model, ch):
    """Forward one-body propagators bu, bd for one recorded step (B = expK D expK)."""
    n = model.n
    du = np.ones(n); dd = np.ones(n)
    for (a, sa, b, sb, e1, e2, cf), k in zip(model.terms, ch):
        (du if sa == 0 else dd)[a] *= e1[k]
        (du if sb == 0 else dd)[b] *= e2[k]
    bu = model.expK @ (du[:, None] * model.expK)
    bd = model.expK_dn @ (dd[:, None] * model.expK_dn)
    return bu, bd


def measure_bp(model, ket_up, ket_dn, rec, w, Phi_u, Phi_d, Fd, etas):
    """BP-AGP d-wave vertex + energy at each eta, from the recorded window. Returns
    per-eta accumulators (weighted sums), total SIGNED weight, and total |weight|
    (for <sign> under constraint release)."""
    n = model.n; I = np.eye(n)
    acc = {e: dict(E=0.0, Sfull=0.0, Pu=np.zeros((n, n)), Pd=np.zeros((n, n)),
                   FKd=0.0, FKa=0.0) for e in etas}
    sw = 0.0; saw = 0.0
    for i in range(len(w)):
        if w[i] == 0 or any(c is None for c in rec[i]):
            continue
        # forward propagators over the window:  Bu_tot = bu_L ... bu_1
        Bu = np.eye(n); Bd = np.eye(n)
        for ch in rec[i]:
            bu, bd = _step_bmats(model, ch)
            Bu = bu @ Bu; Bd = bd @ Bd
        A0 = ket_up[i]; B0 = ket_dn[i]
        for e in etas:
            F = augmented_geminal(Phi_u, Phi_d, Fd, e)
            F_bp = Bu.T @ F @ Bd                       # BP-AGP geminal seen at slice 0
            try:
                c = agp_contractions(F_bp, A0, B0)
            except np.linalg.LinAlgError:
                continue
            Gu, Gd, Kd, Ka = c["Gu"], c["Gd"], c["Kd"], c["Ka"]
            Pu = I - Gu.T; Pd = I - Gd.T
            normal = float((Pu * (Fd @ Pd @ Fd.T)).sum())
            FKd = float((Fd * Kd).sum()); FKa = float((Fd * Ka).sum())
            a = acc[e]
            a["E"] += w[i] * agp_local_energy(model, F_bp, A0, B0)
            a["Sfull"] += w[i] * (normal - FKa * FKd)
            a["Pu"] += w[i] * Pu; a["Pd"] += w[i] * Pd
            a["FKd"] += w[i] * FKd; a["FKa"] += w[i] * FKa
        sw += w[i]; saw += abs(w[i])
    return acc, sw, saw


def comb(model, pu, pd, w, rng):
    nw = len(w)
    for i in range(nw):
        if w[i] > 0:
            pu[i] = np.linalg.qr(pu[i])[0]; pd[i] = np.linalg.qr(pd[i])[0]
    wc = np.clip(w, 0, None); s = wc.sum()
    if s > 0:
        wc *= nw / s; cum = np.cumsum(wc) / wc.sum()
        idx = np.clip(np.searchsorted(cum, (rng.random() / nw) + np.arange(nw) / nw), 0, nw - 1)
        pu = pu[idx].copy(); pd = pd[idx].copy(); w = np.ones(nw)
    return pu, pd, w


def one_pop(P, seed):
    Ku, Kd = am_hopping(P["lx"], P["ly"], P["t0"], P["tam"], P["t1"])
    model = LatticeModel(Ku, P["lx"], P["ly"], P["dt"], P["U"], 0.0, 0.0, K_dn=Kd)
    vu, vd = model.eigvecs, model.eigvecs_dn
    Phi_u = vu[:, :P["nup"]].copy(); Phi_d = vd[:, :P["ndn"]].copy()
    _, Fd = build_Fd(P["lx"], P["ly"])
    etas = P["etas"]; bp = P["bp"]; nw = P["nw"]
    rng = np.random.default_rng(seed)
    pu = np.stack([Phi_u.copy() for _ in range(nw)])
    pd = np.stack([Phi_d.copy() for _ in range(nw)])
    w = np.ones(nw)
    # equilibrate (free constraint, no recording)
    for it in range(P["nequil"]):
        rec = [[] for _ in range(nw)]
        pu, pd, w = step_record_free(model, Phi_u, Phi_d, pu, pd, w, rng, rec)
        if (it + 1) % 10 == 0:
            pu, pd, w = comb(model, pu, pd, w, rng)
    # measurement blocks: ket at slice 0, record bp forward steps, BP-AGP measure.
    # release=True: the measurement window uses FREE (signed) projection on a COPY of
    # the CP ensemble (removing the constraint/ket bias), while the MAIN walkers keep
    # evolving under CP for decorrelation. release=False: the window is CP (pure BP).
    rel = P["release"]
    tot = {e: dict(E=0.0, Sfull=0.0, Pu=np.zeros((model.n, model.n)),
                   Pd=np.zeros((model.n, model.n)), FKd=0.0, FKa=0.0) for e in etas}
    W = 0.0; SW = 0.0; SAW = 0.0
    for blk in range(P["nblocks"]):
        mpu = pu.copy(); mpd = pd.copy(); mw = w.copy()    # release copy (ket at slice 0)
        ket_up = mpu.copy(); ket_dn = mpd.copy()
        rec = [[] for _ in range(nw)]
        for l in range(bp):
            mpu, mpd, mw = step_record_free(model, Phi_u, Phi_d, mpu, mpd, mw, rng, rec,
                                            release=rel)
            if rel:                                        # rescale signed weights (cancels in ratios)
                mag = np.mean(np.abs(mw[mw != 0])) if np.any(mw != 0) else 1.0
                if mag > 0: mw = mw / mag
        acc, sw, saw = measure_bp(model, ket_up, ket_dn, rec, mw, Phi_u, Phi_d, Fd, etas)
        if saw > 0:
            for e in etas:
                for key in ("E", "Sfull", "Pu", "Pd", "FKd", "FKa"):
                    tot[e][key] += acc[e][key]
            W += sw; SW += sw; SAW += saw
        # advance MAIN ensemble under CP (+comb) for the next block
        for l in range(bp):
            rec2 = [[] for _ in range(nw)]
            pu, pd, w = step_record_free(model, Phi_u, Phi_d, pu, pd, w, rng, rec2)
        pu, pd, w = comb(model, pu, pd, w, rng)
    out = {"_sign": (SW / SAW) if SAW > 0 else np.nan}
    for e in etas:
        a = tot[e]
        if abs(W) < 1e-300:
            out[e] = dict(E=np.nan, Sfull=np.nan, vtx=np.nan); continue
        E = a["E"] / W; Sfull = a["Sfull"] / W
        Pu = a["Pu"] / W; Pd = a["Pd"] / W
        bubble = float((Pu * (Fd @ Pd @ Fd.T)).sum())
        # connected vertex subtracts BOTH disconnected pieces: the normal bubble AND
        # the anomalous pair-amplitude product <Delta><Delta^dag> (nonzero for the
        # paired AGP bra). vtx = <S_full> - bubble_normal + <Fd:Ka><Fd:Kd>. At U=0
        # (no walker fluctuation) this is exactly 0, as it must be.
        anom_disc = (a["FKa"] / W) * (a["FKd"] / W)
        out[e] = dict(E=E, Sfull=Sfull, vtx=Sfull - bubble + anom_disc)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=2); ap.add_argument("--ly", type=int, default=2)
    ap.add_argument("--nup", type=int, default=2); ap.add_argument("--ndn", type=int, default=2)
    ap.add_argument("--U", type=float, default=4.0); ap.add_argument("--t0", type=float, default=1.0)
    ap.add_argument("--tam", type=float, default=0.0); ap.add_argument("--t1", type=float, default=0.3)
    ap.add_argument("--etas", type=str, default="0.0,0.3,0.6,1.0,2.0")
    ap.add_argument("--dt", type=float, default=0.02); ap.add_argument("--nw", type=int, default=60)
    ap.add_argument("--nequil", type=int, default=120); ap.add_argument("--bp", type=int, default=20)
    ap.add_argument("--nblocks", type=int, default=20); ap.add_argument("--npop", type=int, default=16)
    ap.add_argument("--seed", type=int, default=1); ap.add_argument("--ed", action="store_true")
    ap.add_argument("--release", action="store_true",
                    help="free (signed) projection in the measurement window -> removes "
                         "the constraint/ket bias (decaying sign); else pure CP back-prop.")
    ap.add_argument("-o", "--out", type=str, default=None,
                    help="write a JSON record (params + per-eta energy/vertex + ED) for the platform.")
    a = ap.parse_args()
    assert a.nup == a.ndn
    etas = [float(x) for x in a.etas.split(",")]
    P = dict(lx=a.lx, ly=a.ly, nup=a.nup, ndn=a.ndn, U=a.U, t0=a.t0, tam=a.tam, t1=a.t1,
             etas=etas, dt=a.dt, nw=a.nw, nequil=a.nequil, bp=a.bp, nblocks=a.nblocks,
             release=a.release)
    import multiprocessing as mp
    with mp.Pool(min(a.npop, mp.cpu_count())) as pool:
        res = pool.starmap(one_pop, [(P, a.seed + p) for p in range(a.npop)])
    print(f"# {a.lx}x{a.ly} nup={a.nup} ndn={a.ndn} U={a.U} tam={a.tam} t1={a.t1} "
          f"dt={a.dt} bp={a.bp} nblocks={a.nblocks} npop={a.npop} nw={a.nw}")
    mode = "RELEASE window" if a.release else "CP back-prop"
    sgn = np.nanmean([r["_sign"] for r in res])
    print(f"# free-electron constraint + AGP bra; measurement = {mode}; <sign>={sgn:.4f}")
    Eref = Sf = Sv = None
    if a.ed:
        Eref, Sf, Sv = ed_reference(P)
    print(f"{'eta':>6} {'energy_bp':>11} {'S_d full':>10} {'S_d VERTEX':>12} {'+/-':>8}")
    rows = []
    for e in etas:
        E = np.array([r[e]["E"] for r in res]); V = np.array([r[e]["vtx"] for r in res])
        Sfu = np.array([r[e]["Sfull"] for r in res])
        nE = np.sum(np.isfinite(E)); nV = np.sum(np.isfinite(V))
        Em = float(np.nanmean(E)); Ee = float(np.nanstd(E) / np.sqrt(max(nE, 1)))
        m = float(np.nanmean(V)); err = float(np.nanstd(V) / np.sqrt(max(nV, 1)))
        print(f"{e:>6.2f} {Em:>11.4f} {np.nanmean(Sfu):>10.3f} {m:>12.4f} {err:>8.4f}")
        rows.append(dict(eta=e, energy=Em, energy_err=Ee, S_full=float(np.nanmean(Sfu)),
                         vertex=m, vertex_err=err))
    if a.ed:
        print(f"\n# full-Fock ED: energy={Eref:+.5f}  S_d full={Sf:+.4f}  S_d VERTEX={Sv:+.4f}")
    if a.out:
        import json
        rec = dict(lx=a.lx, ly=a.ly, nup=a.nup, ndn=a.ndn, U=a.U, t0=a.t0, tam=a.tam,
                   t1=a.t1, dt=a.dt, bp=a.bp, nblocks=a.nblocks, npop=a.npop, nw=a.nw,
                   sign=float(sgn), rows=rows,
                   ed=(dict(energy=float(Eref), S_full=float(Sf), S_vertex=float(Sv))
                       if a.ed else None))
        with open(a.out, "w") as fh:
            json.dump(rec, fh, indent=2)


if __name__ == "__main__":
    main()
