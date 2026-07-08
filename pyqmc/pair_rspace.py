#!/usr/bin/env python3
"""Real-space-resolved equal-time d-wave pairing VERTEX N_pp(R), matching the
manuscript's plot_appeal.ipynb (Vertex_Nppair.dat -> N_pp(R=0) and <N_pp(R)>_{|R|>2}).

The manuscript reports the d_{x^2-y^2} connected vertex resolved by the separation R
between two pairs, and looks at the LOCAL (R=0) and LONG-RANGE (|R|>2, the SC signal)
parts -- NOT the q=0 total sum we computed before. This driver reproduces that
observable in pyqmc (free-electron trial, back-propagated), so we can compare to the
Fortran on a matched footing: U=4, closed-shell filling near half (no open-shell
degeneracy), tam scan.

Per (m,n): connected pair correlation
    P_d(m,n) = <Delta_d^dag(m) Delta_d(n)>_conn,  Delta_d^dag(m)=sum_delta f_d c^+_{m up} c^+_{m+delta dn}
             = Gu[m,n] (Fd Gd Fd^T)[m,n]  -  bubble(m,n)   (bubble from ensemble-avg GFs)
binned by R = r_m - r_n (min-image PBC). Reports N_pp(R=0), <N_pp>_{|R|>2}, and the
q=0 sum (= old observable) for cross-reference.

    source tools/env.sh
    python pyqmc/pair_rspace.py --lx 8 --ly 8 --nup 31 --ndn 31 --U 4 --tam 0.3   # closed shell
    # U=0 gate: vertex ~ 0 at all R
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
import argparse, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from cpqmc import LatticeModel, am_hopping, _green
from agp_dwave import build_Fd
from agp_bp_vertex import step_record_free, comb


def bp_bra_free(model, Lu0, Ld0, rec_i, bp):
    """Back-propagate the free trial through the recorded window (mirror
    Estimators.bp_bra): apply e^{-dtK/2} D e^{-dtK/2} in reverse time order."""
    Lu = Lu0.copy(); Ld = Ld0.copy()
    for s in range(bp - 1, -1, -1):
        ch = rec_i[s]
        if ch is None:
            return None, None
        Lu = model.expK @ Lu; Ld = model.expK_dn @ Ld
        for (a, sa, b, sb, e1, e2, cf), k in zip(model.terms, ch):
            if sa == 0: Lu[a, :] *= e1[k]
            else:       Ld[a, :] *= e1[k]
            if sb == 0: Lu[b, :] *= e2[k]
            else:       Ld[b, :] *= e2[k]
        Lu = model.expK @ Lu; Ld = model.expK_dn @ Ld
    return Lu, Ld


def sep_index(lx, ly):
    """For each (m,n) site pair return the min-image |R| (Euclidean) and a key for
    R=0. Sites indexed i = x*ly + y (am_hopping convention)."""
    n = lx * ly
    xs = np.arange(n) // ly; ys = np.arange(n) % ly
    dx = (xs[:, None] - xs[None, :]) % lx; dx = np.minimum(dx, lx - dx)
    dy = (ys[:, None] - ys[None, :]) % ly; dy = np.minimum(dy, ly - dy)
    R = np.sqrt(dx ** 2 + dy ** 2)
    return R


def one_pop(P, seed):
    Ku, Kd = am_hopping(P["lx"], P["ly"], P["t0"], P["tam"], P["t1"])
    model = LatticeModel(Ku, P["lx"], P["ly"], P["dt"], P["U"], 0.0, 0.0, K_dn=Kd)
    n = model.n
    Phi_u = model.eigvecs[:, :P["nup"]].copy(); Phi_d = model.eigvecs_dn[:, :P["ndn"]].copy()
    _, Fd = build_Fd(P["lx"], P["ly"])
    rng = np.random.default_rng(seed); nw = P["nw"]; bp = P["bp"]
    pu = np.stack([Phi_u.copy() for _ in range(nw)]); pd = np.stack([Phi_d.copy() for _ in range(nw)])
    w = np.ones(nw)
    for _ in range(P["nequil"]):
        rec0 = [[] for _ in range(nw)]
        pu, pd, w = step_record_free(model, Phi_u, Phi_d, pu, pd, w, rng, rec0)
        pu, pd, w = comb(model, pu, pd, w, rng)
    Pfull = np.zeros((n, n)); Gu_s = np.zeros((n, n)); Gd_s = np.zeros((n, n)); W = 0.0
    for _ in range(P["nblocks"]):
        ket_up = pu.copy(); ket_dn = pd.copy()
        rec = [[] for _ in range(nw)]
        for _ in range(bp):
            pu, pd, w = step_record_free(model, Phi_u, Phi_d, pu, pd, w, rng, rec)
        for i in range(nw):
            if w[i] <= 0 or any(c is None for c in rec[i]):
                continue
            Lu, Ld = bp_bra_free(model, Phi_u, Phi_d, rec[i], bp)
            if Lu is None:
                continue
            try:
                Gu = _green(Lu, ket_up[i]).T; Gd = _green(Ld, ket_dn[i]).T  # [i,j]=<c+_i c_j>
            except np.linalg.LinAlgError:
                continue
            Pfull += w[i] * (Gu * (Fd @ Gd @ Fd.T))
            Gu_s += w[i] * Gu; Gd_s += w[i] * Gd; W += w[i]
        pu, pd, w = comb(model, pu, pd, w, rng)
    if W <= 0:
        return None
    Pfull /= W; Gu_a = Gu_s / W; Gd_a = Gd_s / W
    bubble = Gu_a * (Fd @ Gd_a @ Fd.T)
    Pvtx = Pfull - bubble                      # connected vertex matrix N_pp(m,n)
    return Pvtx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=8); ap.add_argument("--ly", type=int, default=8)
    ap.add_argument("--nup", type=int, default=31); ap.add_argument("--ndn", type=int, default=31)
    ap.add_argument("--U", type=float, default=4.0); ap.add_argument("--t0", type=float, default=1.0)
    ap.add_argument("--tam", type=float, default=0.0); ap.add_argument("--t1", type=float, default=0.0)
    ap.add_argument("--dt", type=float, default=0.05); ap.add_argument("--nw", type=int, default=120)
    ap.add_argument("--nequil", type=int, default=120); ap.add_argument("--nblocks", type=int, default=40)
    ap.add_argument("--bp", type=int, default=16); ap.add_argument("--npop", type=int, default=12)
    ap.add_argument("--seed", type=int, default=1); ap.add_argument("--rmin", type=float, default=2.0)
    ap.add_argument("--csv", action="store_true", help="one line: lx ly nup ndn U tam t1 seed Npp_R0 Npp_Rgt Npp_qsum")
    a = ap.parse_args()
    P = dict(lx=a.lx, ly=a.ly, nup=a.nup, ndn=a.ndn, U=a.U, t0=a.t0, tam=a.tam, t1=a.t1,
             dt=a.dt, nw=a.nw, nequil=a.nequil, nblocks=a.nblocks, bp=a.bp)
    import multiprocessing as mp
    with mp.Pool(min(a.npop, mp.cpu_count())) as pool:
        mats = [m for m in pool.starmap(one_pop, [(P, a.seed + p) for p in range(a.npop)]) if m is not None]
    R = sep_index(a.lx, a.ly)
    r0 = (R == 0.0); rgt = (R > a.rmin)
    r0v = np.array([m[r0].mean() for m in mats])
    rgtv = np.array([m[rgt].mean() for m in mats])
    qsum = np.array([m.sum() for m in mats])      # q=0 total (old observable)
    def ms(x): return x.mean(), (x.std() / np.sqrt(len(x)) if len(x) > 1 else 0.0)
    if a.csv:
        print("%d %d %d %d %g %g %g %d %.6f %.6f %.6f %.6f %.6f %.6f" % (
            a.lx, a.ly, a.nup, a.ndn, a.U, a.tam, a.t1, a.seed,
            *ms(r0v), *ms(rgtv), *ms(qsum)), flush=True)
        return
    print(f"# {a.lx}x{a.ly} nup={a.nup} ndn={a.ndn} n={2*a.nup/(a.lx*a.ly):.3f} U={a.U} tam={a.tam} t1={a.t1} npop={len(mats)}")
    m0, e0 = ms(r0v); mg, eg = ms(rgtv); mq, eq = ms(qsum)
    print(f"N_pp(R=0)        = {m0:.4f} +/- {e0:.4f}")
    print(f"<N_pp>_(|R|>{a.rmin:g})   = {mg:.5f} +/- {eg:.5f}   <- long-range SC diagnostic")
    print(f"q=0 sum (old)    = {mq:.4f} +/- {eq:.4f}")


if __name__ == "__main__":
    main()
