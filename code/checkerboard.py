#!/usr/bin/env python3
"""CHECKERBOARD (planar-pyrochlore) single-band Hubbard for the qmc-platform
susceptibility engine, WITHOUT touching the validated cpqmc.py.

Adds three things as standalone functions operating on a CPMC/Estimators instance:
  1. checkerboard_hopping(...)  -- the spin-independent K matrix, IDENTICAL to
     Checkerboard_Model/checkerboard_ed.build_hopping == mc2duph.f90/GetK.
  2. the dxy diagonal-bond pair form factor (the channel cpqmc lacks; s + dx2-y2
     are the only built-ins there).
  3. run_bp_chid_cb(...) -- the unequal-time singlet PAIRING susceptibility
     chi_a = int C_a(tau) dtau (full + connected/vertex) for a = {s, dx2-y2, dxy}.

Convention K = hopping matrix in H = sum K c+c (production t0 = -1 -> physical NN
hop +1). t1 = -t' (=+0.3), t2 = -delta. AM is EMERGENT (spin-independent K).
"""
from __future__ import annotations
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from cpqmc import _green   # module-level Green's helper


# ----------------------------------------------------------------------
# 1. checkerboard hopping (matches checkerboard_ed.build_hopping / GetK)
# ----------------------------------------------------------------------
def checkerboard_hopping(lx, ly, t0=-1.0, t1=0.3, t2=0.0, apx=1, apy=1):
    """nsites x nsites checkerboard one-body matrix (spin-independent).
    A sites (x+y even): '/'=t1+t2, '\\'=t1-t2 ; B sites (x+y odd): swapped.
    For delta: pass t2 = -delta. Feed as CPMC(..., K=this, K_dn=None).

    apx, apy in {+1 (periodic, default), -1 (anti-periodic)}: boundary phase in x/y.
    Anti-periodic (-1) stays REAL (so CPQMC works unchanged) and shifts the k-grid off
    (0,pi),(pi,0), opening a gap at HALF-FILLING (open-shell under periodic BC). apx=apy=1
    reproduces the original periodic matrix bit-for-bit."""
    n = lx * ly
    K = np.zeros((n, n))
    tp = t1 + t2
    tm = t1 - t2
    def idx(x, y): return (x % lx) * ly + (y % ly)
    def bs(x, y, dx, dy):                     # boundary sign for a hop crossing the edge
        sx = apx if (x + dx < 0 or x + dx >= lx) else 1
        sy = apy if (y + dy < 0 or y + dy >= ly) else 1
        return sx * sy
    for x in range(lx):
        for y in range(ly):
            i = idx(x, y)
            K[i, idx(x + 1, y)] = t0 * bs(x, y, 1, 0)
            K[i, idx(x - 1, y)] = t0 * bs(x, y, -1, 0)
            K[i, idx(x, y + 1)] = t0 * bs(x, y, 0, 1)
            K[i, idx(x, y - 1)] = t0 * bs(x, y, 0, -1)
            a, b = (tp, tm) if (x + y) % 2 == 0 else (tm, tp)
            K[i, idx(x + 1, y + 1)] = a * bs(x, y, 1, 1)    # '/'  main diagonal
            K[i, idx(x - 1, y - 1)] = a * bs(x, y, -1, -1)
            K[i, idx(x + 1, y - 1)] = b * bs(x, y, 1, -1)   # '\'  anti-diagonal
            K[i, idx(x - 1, y + 1)] = b * bs(x, y, -1, 1)
    return 0.5 * (K + K.T)


# ----------------------------------------------------------------------
# 2. bond form factors (single-band square index i = x*ly + y)
# ----------------------------------------------------------------------
def nn_bond_factors(lx, ly):
    """NN s-wave (all +1) and d_{x^2-y^2} (+1 on x, -1 on y). Copy of
    cpqmc.CPMC._bond_factors for the single band, so results are directly
    comparable to run_bp_chid."""
    n = lx * ly
    Fs = np.zeros((n, n)); Fd = np.zeros((n, n))
    def idx(x, y): return (x % lx) * ly + (y % ly)
    for x in range(lx):
        for y in range(ly):
            m = idx(x, y)
            for (j, fd) in [(idx(x + 1, y), +1.0), (idx(x - 1, y), +1.0),
                            (idx(x, y + 1), -1.0), (idx(x, y - 1), -1.0)]:
                Fs[m, j] += 1.0
                Fd[m, j] += fd
    return Fs, Fd


def diag_bond_factors(lx, ly):
    """dxy pair form factor on the DIAGONAL bonds: +1 on the main diagonal
    (+-(1,1)) and -1 on the anti-diagonal (+-(1,-1)) -> f(k) ~ sin kx sin ky."""
    n = lx * ly
    Fdxy = np.zeros((n, n))
    def idx(x, y): return (x % lx) * ly + (y % ly)
    for x in range(lx):
        for y in range(ly):
            m = idx(x, y)
            for (j, f) in [(idx(x + 1, y + 1), +1.0), (idx(x - 1, y - 1), +1.0),
                           (idx(x - 1, y + 1), -1.0), (idx(x + 1, y - 1), -1.0)]:
                Fdxy[m, j] += f
    return Fdxy


# ----------------------------------------------------------------------
# 3. multi-channel unequal-time pairing susceptibility
# ----------------------------------------------------------------------
def _chid_block_multi(est, walkers, ket_up, ket_dn, rec, bp, Ffac):
    """Generalization of Estimators.chid_block to an arbitrary dict of bond form
    factors Ffac = {tag: F}. Returns (C{tag:(L,)}, Pu_s(L,n,n), Pd_s(L,n,n), W).
    C is the FULL (bubble+vertex) per-walker sum; the caller builds the bubble
    from the ensemble-averaged Pu_s/Pd_s exactly as run_bp_chid does."""
    L = bp + 1
    n = est.m.n
    C = {tag: np.zeros(L) for tag in Ffac}
    W = 0.0
    Pu_s = np.zeros((L, n, n)); Pd_s = np.zeros((L, n, n))
    for i in range(walkers.nw):
        if walkers.w[i] <= 0:
            continue
        chs = rec[i]
        if any(c is None for c in chs):
            continue
        Lu, Ld = est.bp_bra(chs, bp)
        if Lu is None:
            continue
        try:
            gu = _green(Lu, ket_up[i]).T
            gd = _green(Ld, ket_dn[i]).T
        except np.linalg.LinAlgError:
            continue
        w = walkers.w[i]
        Bu = np.eye(n); Bd = np.eye(n)
        ImguT = np.eye(n) - gu.T; ImgdT = np.eye(n) - gd.T
        for l in range(L):
            if l > 0:
                bu, bd, _, _ = est._step_bmats(chs[l - 1])
                Bu = bu @ Bu; Bd = bd @ Bd
            Pu = Bu @ ImguT; Pd = Bd @ ImgdT           # time-displaced particle GFs
            for tag, F in Ffac.items():
                C[tag][l] += w * float(np.sum(Pu * (F @ Pd @ F.T)))
            Pu_s[l] += w * Pu; Pd_s[l] += w * Pd
        W += w
    return C, Pu_s, Pd_s, W


def run_bp_chid_cb(q, Ffac, nequil=150, nblocks=30, bp=20, ortho=10, pc=10):
    """Multi-channel unequal-time pairing susceptibility on a CPMC instance q.
    chi_a = int C_a(tau) dtau (trapezoidal), full and connected/VERTEX, per
    channel tag in Ffac. Mirrors CPMC.run_bp_chid exactly (block loop, bubble
    from averaged GFs, per-block error bars)."""
    for it in range(nequil):
        q.step()
        if (it + 1) % ortho == 0: q.reorthogonalize()
        if (it + 1) % pc == 0: q.pop_control()
        q._maybe_update_trial(it)
    acc = {f"{tag}::{k}": [] for tag in Ffac for k in ("full", "vertex")}
    for blk in range(nblocks):
        q.reorthogonalize()
        ket_up = q.walkers.phi_up.copy(); ket_dn = q.walkers.phi_dn.copy()
        rec = [[] for _ in range(q.nw)]
        for _ in range(bp):
            q.prop.step_record(q.walkers, rec)
        C, Pu_s, Pd_s, W = _chid_block_multi(q.est, q.walkers, ket_up, ket_dn, rec, bp, Ffac)
        if W <= 0:
            q.pop_control(); continue
        Pu = Pu_s / W; Pd = Pd_s / W                    # ensemble-averaged GFs -> bubble
        for tag, F in Ffac.items():
            Cfull = C[tag] / W
            bubble = np.array([float(np.sum(Pu[l] * (F @ Pd[l] @ F.T)))
                               for l in range(bp + 1)])
            acc[f"{tag}::full"].append(Cfull)
            acc[f"{tag}::vertex"].append(Cfull - bubble)
        q.pop_control()
    taus = q.dt * np.arange(bp + 1)
    out = {"nsites": q.n, "bp": bp, "dt": q.dt, "taus": taus.tolist()}
    for key, blocks in acc.items():
        Bm = np.array(blocks)
        chi_blk = q.dt * (Bm[:, 1:-1].sum(axis=1) + 0.5 * (Bm[:, 0] + Bm[:, -1]))
        tag, kind = key.split("::")
        suffix = "" if kind == "full" else "_vertex"
        out[f"C_{tag}_tau{suffix}"] = Bm.mean(axis=0).tolist()
        out[f"chi_{tag}{suffix}"] = float(chi_blk.mean())
        out[f"chi_{tag}{suffix}_err"] = float(chi_blk.std() / np.sqrt(len(chi_blk))
                                              if len(chi_blk) > 1 else 0.0)
    return out
