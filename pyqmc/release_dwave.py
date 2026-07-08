#!/usr/bin/env python3
"""Constraint RELEASE for the d-wave VERTEX (attack the d-wave CP bias).

Equilibrate under the constrained path (positive weights, no node crossing), then
RELEASE the constraint (free projection, signed weights) and measure the equal-time
d-wave vertex at each release step tau = l*dt. Under release the node bias is
removed and the estimate should climb from the (biased) CP value toward the EXACT
value -- until the average sign <s> decays and the noise blows up. We watch both.

Validation target: 4x2, tam=0, t1=0.3, U=4, half filling -> ED d-wave vertex = +23.18
(CP-AFQMC gives ~ -2). Does release recover the ED value before the sign dies?

    source tools/env.sh
    python pyqmc/release_dwave.py --lx 4 --ly 2 --nup 4 --ndn 4 --U 4 --t1 0.3
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import argparse, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from cpqmc import LatticeModel, TrialWF, am_hopping, _green


def build_F(lx, ly):
    n = lx * ly
    Fs = np.zeros((n, n)); Fd = np.zeros((n, n))
    def idx(x, y): return (x % lx) * ly + (y % ly)
    for x in range(lx):
        for y in range(ly):
            m = idx(x, y)
            for (dx, dy, fd) in ((1, 0, 1.0), (-1, 0, 1.0), (0, 1, -1.0), (0, -1, -1.0)):
                j = idx(x + dx, y + dy); Fs[m, j] += 1.0; Fd[m, j] += fd
    return Fs, Fd


def step_pop(model, trial, pu, pd, w, rng, constrained):
    """One propagation step over the population. constrained: q=max(R,0), positive
    weight (CP). else: q=|R|, signed weight (free projection)."""
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


def measure_vertex(model, trial, pu, pd, w, Fs, Fd):
    n = model.n; sw = 0.0; sa = 0.0; full_d = 0.0; full_s = 0.0
    Gu_a = np.zeros((n, n)); Gd_a = np.zeros((n, n))
    for i in range(len(w)):
        if w[i] == 0:
            continue
        try:
            Gu = _green(trial.up, pu[i]).T; Gd = _green(trial.dn, pd[i]).T
        except np.linalg.LinAlgError:
            continue
        full_d += w[i] * float((Gu * (Fd @ Gd @ Fd.T)).sum())
        full_s += w[i] * float((Gu * (Fs @ Gd @ Fs.T)).sum())
        Gu_a += w[i] * Gu; Gd_a += w[i] * Gd
        sw += w[i]; sa += abs(w[i])
    if sw == 0:
        return np.nan, np.nan, 0.0
    Gu_a /= sw; Gd_a /= sw
    disc_d = float((Gu_a * (Fd @ Gd_a @ Fd.T)).sum())
    return full_d / sw, full_d / sw - disc_d, sw / sa   # full_d, vertex_d, sign


def one_pop(P, seed):
    Ku, Kd = am_hopping(P["lx"], P["ly"], P["t0"], P["tam"], P["t1"])
    model = LatticeModel(Ku, P["lx"], P["ly"], P["dt"], P["U"], 0.0, 0.0, K_dn=Kd)
    vu, vd = model.eigvecs, model.eigvecs_dn
    trial = TrialWF(vu[:, :P["nup"]].copy(), vd[:, :P["ndn"]].copy(), mode="fixed")
    Fs, Fd = build_F(P["lx"], P["ly"])
    rng = np.random.default_rng(seed)
    nw = P["nw"]
    pu = np.stack([trial.up.copy() for _ in range(nw)])
    pd = np.stack([trial.dn.copy() for _ in range(nw)])
    w = np.ones(nw)
    # CP equilibration (positive weights, periodic reorth + simple pop control)
    for it in range(P["nequil"]):
        pu, pd, w = step_pop(model, trial, pu, pd, w, rng, constrained=True)
        if (it + 1) % 10 == 0:
            for i in range(nw):
                if w[i] > 0:
                    pu[i] = np.linalg.qr(pu[i])[0]; pd[i] = np.linalg.qr(pd[i])[0]
            wc = np.clip(w, 0, None); s = wc.sum()
            if s > 0:                                  # comb resample
                wc *= nw / s; cum = np.cumsum(wc) / wc.sum()
                idx = np.clip(np.searchsorted(cum, (rng.random() / nw) + np.arange(nw) / nw), 0, nw - 1)
                pu = pu[idx].copy(); pd = pd[idx].copy(); w = np.ones(nw)
    # RELEASE: signed free projection, no reorth/comb (short), measure each step
    out = [measure_vertex(model, trial, pu, pd, w, Fs, Fd)]   # tau=0 (CP value)
    for it in range(P["nrelease"]):
        pu, pd, w = step_pop(model, trial, pu, pd, w, rng, constrained=False)
        mag = np.mean(np.abs(w[w != 0])) if np.any(w != 0) else 1.0
        if mag > 0:
            w = w / mag                                # rescale (cancels in ratios)
        out.append(measure_vertex(model, trial, pu, pd, w, Fs, Fd))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=4); ap.add_argument("--ly", type=int, default=2)
    ap.add_argument("--nup", type=int, default=4); ap.add_argument("--ndn", type=int, default=4)
    ap.add_argument("--U", type=float, default=4.0); ap.add_argument("--t0", type=float, default=1.0)
    ap.add_argument("--tam", type=float, default=0.0); ap.add_argument("--t1", type=float, default=0.3)
    ap.add_argument("--dt", type=float, default=0.02); ap.add_argument("--nw", type=int, default=80)
    ap.add_argument("--nequil", type=int, default=120); ap.add_argument("--nrelease", type=int, default=18)
    ap.add_argument("--npop", type=int, default=24); ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    P = dict(lx=a.lx, ly=a.ly, nup=a.nup, ndn=a.ndn, U=a.U, t0=a.t0, tam=a.tam, t1=a.t1,
             dt=a.dt, nw=a.nw, nequil=a.nequil, nrelease=a.nrelease)
    import multiprocessing as mp
    with mp.Pool(min(a.npop, mp.cpu_count())) as pool:
        res = pool.starmap(one_pop, [(P, a.seed + p) for p in range(a.npop)])
    L = a.nrelease + 1
    print(f"# {a.lx}x{a.ly} nup={a.nup} ndn={a.ndn} U={a.U} tam={a.tam} t1={a.t1} dt={a.dt}  npop={a.npop}")
    print(f"# constraint release of the d-wave VERTEX. tau=0 is the CP value.")
    print(f"{'tau':>6} {'<sign>':>8} {'d_full':>9} {'d_VERTEX':>10} {'+/-':>8}")
    for l in range(L):
        vts = np.array([r[l][1] for r in res]); fds = np.array([r[l][0] for r in res])
        sgn = np.array([r[l][2] for r in res])
        vts = vts[np.isfinite(vts)]
        m = vts.mean(); e = vts.std() / np.sqrt(len(vts)) if len(vts) > 1 else np.nan
        print(f"{l * a.dt:>6.3f} {np.nanmean(sgn):>8.4f} {np.nanmean(fds):>9.3f} {m:>10.3f} {e:>8.3f}")


if __name__ == "__main__":
    main()
