#!/usr/bin/env python3
"""UNIFIED single-walk d-wave estimator with FULL k-resolution. One CPMC back-prop
pass measures, for the d_{x^2-y^2} channel, BOTH the equal-time pair CORRELATION and
the tau-integrated SUSCEPTIBILITY, each in FULL and connected-VERTEX form, fully
resolved in pair center-of-mass momentum q (and real-space separation R). From the
R-vector-resolved S(R) [= sum_m <Delta^dag(m) Delta(m+R)>] we report four reductions:

  maxk : max over q of P(q)   (peak pairing momentum)
  k0   : P(q=0) = sum_R S(R)  (uniform channel -- old 'q=0 sum')
  r0   : S(R=0)/N             (local/onsite)
  rgt  : <S(R)/N>_{|R|>2}     (long-range, the SC diagnostic)

So one parameter point yields correlation & susceptibility x {full,vtx} x {maxk,k0,r0,rgt}.
Model: am_hopping(t0,tam,t1,tp); t1 = spin-DEPENDENT (anisotropic, altermagnet) NNN,
tp = spin-INDEPENDENT (isotropic, standard) NNN -- so we can separate the two.

    python pyqmc/unified_scan.py --lx 8 --ly 8 --nup 31 --ndn 31 --U 4 --tam 0.3
    python pyqmc/unified_scan.py ... --t1 0.2     # anisotropic NNN
    python pyqmc/unified_scan.py ... --tp 0.2     # isotropic  NNN
    python pyqmc/unified_scan.py ... --csv        # machine-readable (col order in code)
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
from pair_rspace import bp_bra_free


def step_bmats(model, ch):
    n = model.n; d_up = np.ones(n); d_dn = np.ones(n)
    for (a, sa, b, sb, e1, e2, cf), k in zip(model.terms, ch):
        (d_up if sa == 0 else d_dn)[a] *= e1[k]
        (d_up if sb == 0 else d_dn)[b] *= e2[k]
    bu = model.expK @ (d_up[:, None] * model.expK)
    bd = model.expK_dn @ (d_dn[:, None] * model.expK_dn)
    return bu, bd


def _shift_index(lx, ly):
    """shift[dx,dy] = array n(m) mapping site m -> site (x+dx, y+dy) (PBC). Site i=x*ly+y."""
    n = lx * ly; xm = np.arange(n) // ly; ym = np.arange(n) % ly
    S = np.empty((lx, ly, n), dtype=int)
    for dx in range(lx):
        for dy in range(ly):
            S[dx, dy] = ((xm + dx) % lx) * ly + (ym + dy) % ly
    return S


def reduce_mat(mat, shift, rmin):
    """matrix M[m,n] = <Delta^dag(m) Delta(n)> -> S(R)=sum_m M[m, m+R] (lx,ly grid),
    then P(q)=FFT2(S). Returns dict(maxk, k0, r0, rgt)."""
    lx, ly, n = shift.shape
    rows = np.arange(n)
    Sg = np.array([[mat[rows, shift[dx, dy]].sum() for dy in range(ly)] for dx in range(lx)])
    Pq = np.real(np.fft.fft2(Sg))
    # |R| grid (min-image)
    dxg = np.minimum(np.arange(lx), lx - np.arange(lx))[:, None]
    dyg = np.minimum(np.arange(ly), ly - np.arange(ly))[None, :]
    Rmag = np.sqrt(dxg ** 2 + dyg ** 2) * np.ones((lx, ly))
    rgt = (Sg[Rmag > rmin] / n).mean() if np.any(Rmag > rmin) else np.nan
    return dict(maxk=float(Pq.max()), k0=float(Sg.sum()), r0=float(Sg[0, 0] / n), rgt=float(rgt))


def one_pop(P, seed):
    Ku, Kd = am_hopping(P["lx"], P["ly"], P["t0"], P["tam"], P["t1"], P["tp"])
    model = LatticeModel(Ku, P["lx"], P["ly"], P["dt"], P["U"], 0.0, 0.0, K_dn=Kd)
    n = model.n; I = np.eye(n); L = P["bp"] + 1; dt = P["dt"]
    tw = np.full(L, dt); tw[0] = tw[-1] = 0.5 * dt          # trapezoid weights
    Phi_u = model.eigvecs[:, :P["nup"]].copy(); Phi_d = model.eigvecs_dn[:, :P["ndn"]].copy()
    Fs, Fd = build_Fd(P["lx"], P["ly"])
    rng = np.random.default_rng(seed); nw = P["nw"]; bp = P["bp"]
    pu = np.stack([Phi_u.copy() for _ in range(nw)]); pd = np.stack([Phi_d.copy() for _ in range(nw)])
    w = np.ones(nw)
    for _ in range(P["nequil"]):
        rec0 = [[] for _ in range(nw)]
        pu, pd, w = step_record_free(model, Phi_u, Phi_d, pu, pd, w, rng, rec0)
        pu, pd, w = comb(model, pu, pd, w, rng)
    corrF = np.zeros((n, n)); suscF = np.zeros((n, n))     # full matrices (d-wave)
    corrF_s = np.zeros((n, n)); suscF_s = np.zeros((n, n)) # ext-s channel
    Gu_a = np.zeros((n, n)); Gd_a = np.zeros((n, n))        # equal-time bubble avg
    Pu_a = np.zeros((L, n, n)); Pd_a = np.zeros((L, n, n))  # tau bubble avg
    Wd = idxm = Pk = None                                   # tau-integrated k-space pairing eigenvalue
    if P.get("paireig"):
        lx, ly = P["lx"], P["ly"]; xs = np.arange(n) // ly; ys = np.arange(n) % ly
        kxf = (2*np.pi*(np.arange(lx))/lx); kyf = (2*np.pi*(np.arange(ly))/ly)
        KX, KY = np.meshgrid(kxf, kyf, indexing="ij"); kxf = KX.ravel(); kyf = KY.ravel()
        Wd = np.exp(-1j*(np.outer(kxf, xs) + np.outer(kyf, ys)))/np.sqrt(n)
        idxm = np.array([((-(i//ly)) % lx)*ly + ((-(i % ly)) % ly) for i in range(n)])
        Pk = np.zeros((n, n), complex)
    Ek = 0.0; W = 0.0
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
                gu = _green(Lu, ket_up[i]).T; gd = _green(Ld, ket_dn[i]).T
            except np.linalg.LinAlgError:
                continue
            wi = w[i]
            corrF += wi * (gu * (Fd @ gd @ Fd.T))
            corrF_s += wi * (gu * (Fs @ gd @ Fs.T))
            Gu_a += wi * gu; Gd_a += wi * gd
            Ek += wi * (np.sum(model.K * gu) + np.sum(model.K_dn * gd)
                        + P["U"] * np.sum(np.diag(gu) * np.diag(gd)))
            ImguT = I - gu.T; ImgdT = I - gd.T
            Bu = np.eye(n); Bd = np.eye(n)
            for l in range(L):
                if l > 0:
                    bu, bd = step_bmats(model, rec[i][l - 1]); Bu = bu @ Bu; Bd = bd @ Bd
                Pu = Bu @ ImguT; Pd = Bd @ ImgdT
                suscF += (wi * tw[l]) * (Pu * (Fd @ Pd @ Fd.T))
                suscF_s += (wi * tw[l]) * (Pu * (Fs @ Pd @ Fs.T))
                Pu_a[l] += wi * Pu; Pd_a[l] += wi * Pd
                if Wd is not None:                                      # tau-integrated k-pairing matrix
                    Puk = Wd @ Pu @ Wd.conj().T; Pdk = Wd @ Pd @ Wd.conj().T
                    Pk += (wi * tw[l]) * (Puk * Pdk[np.ix_(idxm, idxm)])
            W += wi
        pu, pd, w = comb(model, pu, pd, w, rng)
    if W <= 0:
        return None
    corrF /= W; suscF /= W; corrF_s /= W; suscF_s /= W
    gu_a = Gu_a / W; gd_a = Gd_a / W; Pu_a /= W; Pd_a /= W
    corr_bub = gu_a * (Fd @ gd_a @ Fd.T)
    susc_bub = sum(tw[l] * (Pu_a[l] * (Fd @ Pd_a[l] @ Fd.T)) for l in range(L))
    corr_bub_s = gu_a * (Fs @ gd_a @ Fs.T)
    susc_bub_s = sum(tw[l] * (Pu_a[l] * (Fs @ Pd_a[l] @ Fs.T)) for l in range(L))
    out = dict(energy=Ek / W, corrF=corrF, corrV=corrF - corr_bub,
               suscF=suscF, suscV=suscF - susc_bub,
               corrF_s=corrF_s, corrV_s=corrF_s - corr_bub_s,
               suscF_s=suscF_s, suscV_s=suscF_s - susc_bub_s)
    if Wd is not None:                                     # tau-integrated pairing eigenvalue (T=0, bp window)
        Pk /= W
        Pkbub = sum(tw[l] * ((Wd @ Pu_a[l] @ Wd.conj().T) *
                             (Wd @ Pd_a[l] @ Wd.conj().T)[np.ix_(idxm, idxm)]) for l in range(L))
        Pc = Pk - Pkbub; Pc = 0.5 * (Pc + Pc.conj().T)
        ww, vv = np.linalg.eigh(Pc); phi = vv[:, -1]
        fd = np.cos(kxf) - np.cos(kyf); fs = np.cos(kxf) + np.cos(kyf)
        out["pe_lam"] = float(ww[-1]); out["pe_lam2"] = float(ww[-2])
        out["pe_dov"] = float(abs(np.vdot(fd/np.linalg.norm(fd), phi)))
        out["pe_sov"] = float(abs(np.vdot(fs/np.linalg.norm(fs), phi)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=8); ap.add_argument("--ly", type=int, default=8)
    ap.add_argument("--nup", type=int, default=31); ap.add_argument("--ndn", type=int, default=31)
    ap.add_argument("--U", type=float, default=4.0); ap.add_argument("--t0", type=float, default=1.0)
    ap.add_argument("--tam", type=float, default=0.0); ap.add_argument("--t1", type=float, default=0.0)
    ap.add_argument("--tp", type=float, default=0.0)
    ap.add_argument("--dt", type=float, default=0.05); ap.add_argument("--nw", type=int, default=90)
    ap.add_argument("--nequil", type=int, default=200); ap.add_argument("--nblocks", type=int, default=40)
    ap.add_argument("--bp", type=int, default=16); ap.add_argument("--npop", type=int, default=6)
    ap.add_argument("--seed", type=int, default=1); ap.add_argument("--rmin", type=float, default=2.0)
    ap.add_argument("--csv", action="store_true")
    a = ap.parse_args()
    P = dict(lx=a.lx, ly=a.ly, nup=a.nup, ndn=a.ndn, U=a.U, t0=a.t0, tam=a.tam, t1=a.t1, tp=a.tp,
             dt=a.dt, nw=a.nw, nequil=a.nequil, nblocks=a.nblocks, bp=a.bp)
    import multiprocessing as mp
    with mp.Pool(min(a.npop, mp.cpu_count())) as pool:
        res = [r for r in pool.starmap(one_pop, [(P, a.seed + p) for p in range(a.npop)]) if r]
    shift = _shift_index(a.lx, a.ly)
    fields = ["corrF", "corrV", "suscF", "suscV"]; reds = ["maxk", "k0", "r0", "rgt"]
    out = {}
    for fld in fields:
        for rd in reds:
            vals = np.array([reduce_mat(r[fld], shift, a.rmin)[rd] for r in res])
            out[f"{fld}_{rd}"] = (vals.mean(), vals.std() / np.sqrt(len(vals)) if len(vals) > 1 else 0.0)
    en = np.array([r["energy"] for r in res]); out["energy"] = (en.mean(), en.std() / np.sqrt(len(en)) if len(en) > 1 else 0.0)
    if a.csv:
        order = ["energy"] + [f"{f}_{r}" for f in fields for r in reds]
        flat = " ".join(f"{out[k][0]:.6f}" for k in order)
        print(f"{a.lx} {a.ly} {a.nup} {a.ndn} {a.U:g} {a.tam:g} {a.t1:g} {a.tp:g} {a.seed} {flat}", flush=True)
        return
    print(f"# {a.lx}x{a.ly} N={a.nup}+{a.ndn} n={2*a.nup/(a.lx*a.ly):.3f} U={a.U} "
          f"tam={a.tam} t1={a.t1} tp={a.tp} npop={len(res)}")
    print(f"  energy = {out['energy'][0]:.5f}")
    print(f"  {'':14}{'maxk':>10}{'k0':>10}{'r0':>10}{'rgt':>10}")
    name = {"corrF": "corr full", "corrV": "corr vtx", "suscF": "susc full", "suscV": "susc vtx"}
    for fld in fields:
        print(f"  {name[fld]:<14}" + "".join(f"{out[f'{fld}_{r}'][0]:>10.4f}" for r in reds))


if __name__ == "__main__":
    main()
