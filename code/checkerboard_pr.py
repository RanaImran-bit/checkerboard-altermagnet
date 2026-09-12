#!/usr/bin/env python3
"""Real-space pairing correlations P(R) and their vertex contribution V(R).

This is the Huang, Lin and Gubernatis analysis (PRB 64, 205101): resolve the pair
correlation by the distance R between the two pairs, then split the vertex into a
LOCAL part V(R=0) and a LONG-RANGE part averaged over R > 2. They find extended s
dominates the local piece while d dominates the long-range piece, which is why our
q = 0 sum -- which adds both together -- puts extended s ahead of d_x2-y2 and so
appears to disagree with White.

The existing driver computes

    chi_a = sum over ALL (i,j) of  Pu[i,j] * (F Pd F^T)[i,j]

and discards the site structure. Keeping that matrix and binning entry (i,j) by the
separation between sites i and j gives P(R) instead. Summing P(R) over R must
return the original chi exactly -- checked below and used as the validation.

Distances use the minimum image convention on the periodic L x L lattice, and the
site ordering matches checkerboard_hopping: site s sits at (x, y) = (s // ly, s % ly).

  SIZES=12 US=4 NUPS=72 DELTAS=0,0.4 NPROC=32 NSEED=6 python checkerboard_pr.py
Writes pr_L{L}_U{U}_d{delta}_nup{nup}.csv -- one row per (seed, channel, R).
"""
import os, sys, time, csv
os.environ["OMP_NUM_THREADS"] = "1"; os.environ["MKL_NUM_THREADS"] = "1"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from multiprocessing import Pool
from cpqmc import CPMC
import checkerboard as cb
from checkerboard import _green, _chid_block_multi

T0, T1 = -1.0, 0.3
NW, NEQ, NBLK, BP, DT = 160, 60, int(os.environ.get("NBLK", 40)), 16, 0.05
NPROC = int(os.environ.get("NPROC", min(os.cpu_count(), 32)))
NSEED = int(os.environ.get("NSEED", 6)); SEEDS = list(range(1, NSEED + 1))
SIZES = [int(x) for x in os.environ.get("SIZES", "12").split(",")]
US = [float(x) for x in os.environ.get("US", "4").split(",")]
DELTAS = [float(x) for x in os.environ.get("DELTAS", "0,0.4").split(",")]
NUPS_ENV = os.environ.get("NUPS", "")
N_TARGETS = [2 * k / 36 for k in [9, 10, 12, 14, 16, 18]]
CHAN = ["son", "sext", "d", "dxy"]


def distance_map(lx, ly):
    """R[i,j] = minimum-image distance between sites i and j, and the sorted
    list of distinct distances. Site s is at (s // ly, s % ly)."""
    n = lx * ly
    xs, ys = np.arange(n) // ly, np.arange(n) % ly
    dx = np.abs(xs[:, None] - xs[None, :]); dx = np.minimum(dx, lx - dx)
    dy = np.abs(ys[:, None] - ys[None, :]); dy = np.minimum(dy, ly - dy)
    R = np.sqrt(dx**2 + dy**2)
    Rr = np.round(R, 6)
    return Rr, np.unique(Rr)


def run_point(args):
    """One (L, nup, delta, U, seed). Returns rows of (channel, R, P, bubble, V)."""
    L, nup, delta, U, seed = args
    n = L * L
    K = cb.checkerboard_hopping(L, L, T0, T1, -delta)
    Fs, Fd = cb.nn_bond_factors(L, L); Fdxy = cb.diag_bond_factors(L, L)
    Ffac = {"son": np.eye(n), "sext": Fs, "d": Fd, "dxy": Fdxy}
    q = CPMC(L, L, nup, nup, U=U, dt=DT, nwalkers=NW, seed=seed, K=K, K_dn=None)
    Rmap, Rvals = distance_map(L, L)
    masks = [(Rmap == r) for r in Rvals]

    # Equilibrate exactly as run_bp_chid_cb does. There is no q.run_equil -- an
    # earlier version guarded on hasattr and therefore skipped equilibration
    # entirely, which broke the sum rule by up to 400%.
    ORTHO, PC = 10, 10
    for it in range(NEQ):
        q.step()
        if (it + 1) % ORTHO == 0: q.reorthogonalize()
        if (it + 1) % PC == 0: q.pop_control()
        q._maybe_update_trial(it)

    accP = {c: np.zeros((len(Rvals),)) for c in CHAN}
    accB = {c: np.zeros((len(Rvals),)) for c in CHAN}
    nblk = 0
    for _ in range(NBLK):
        q.reorthogonalize()
        ket_up = q.walkers.phi_up.copy(); ket_dn = q.walkers.phi_dn.copy()
        rec = [[] for _ in range(q.nw)]
        for _ in range(BP):
            q.prop.step_record(q.walkers, rec)
        # per-walker matrices, kept resolved in (i,j) instead of summed
        Msum = {c: np.zeros((BP + 1, n, n)) for c in CHAN}
        Pu_s = np.zeros((BP + 1, n, n)); Pd_s = np.zeros((BP + 1, n, n)); W = 0.0
        for iw in range(q.walkers.nw):
            if q.walkers.w[iw] <= 0:
                continue
            chs = rec[iw]
            if any(c is None for c in chs):
                continue
            Lu, Ld = q.est.bp_bra(chs, BP)
            if Lu is None:
                continue
            try:
                gu = _green(Lu, ket_up[iw]).T
                gd = _green(Ld, ket_dn[iw]).T
            except np.linalg.LinAlgError:
                continue
            w = q.walkers.w[iw]
            Bu = np.eye(n); Bd = np.eye(n)
            ImguT = np.eye(n) - gu.T; ImgdT = np.eye(n) - gd.T
            for l in range(BP + 1):
                if l > 0:
                    bu, bd, _, _ = q.est._step_bmats(chs[l - 1])
                    Bu = bu @ Bu; Bd = bd @ Bd
                Pu = Bu @ ImguT; Pd = Bd @ ImgdT
                for c in CHAN:
                    F = Ffac[c]
                    Msum[c][l] += w * (Pu * (F @ Pd @ F.T))
                Pu_s[l] += w * Pu; Pd_s[l] += w * Pd
            W += w
        if W <= 0:
            q.pop_control(); continue
        Pu_a = Pu_s / W; Pd_a = Pd_s / W
        taus_w = np.full(BP + 1, DT); taus_w[0] = taus_w[-1] = DT / 2   # trapezoid
        for c in CHAN:
            F = Ffac[c]
            Mfull = Msum[c] / W
            Mbub = np.array([Pu_a[l] * (F @ Pd_a[l] @ F.T) for l in range(BP + 1)])
            # integrate over imaginary time first, then bin by separation
            If = np.tensordot(taus_w, Mfull, axes=(0, 0))
            Ib = np.tensordot(taus_w, Mbub, axes=(0, 0))
            accP[c] += np.array([If[m].sum() for m in masks])
            accB[c] += np.array([Ib[m].sum() for m in masks])
        nblk += 1
        q.pop_control()

    rows = []
    for c in CHAN:
        P = accP[c] / max(nblk, 1); B = accB[c] / max(nblk, 1)
        for k, r in enumerate(Rvals):
            rows.append((L, nup, delta, U, seed, c, float(r),
                         float(P[k]), float(B[k]), float(P[k] - B[k])))
    return rows


COLS = ["L", "nup", "delta", "U", "seed", "channel", "R", "P", "bubble", "V"]

if __name__ == "__main__":
    print(f"P(R): sizes={SIZES} US={US} deltas={DELTAS} NSEED={NSEED} on {NPROC} cores",
          flush=True)
    for L in SIZES:
        nups = ([int(x) for x in NUPS_ENV.split(",")] if NUPS_ENV
                else sorted({max(1, round(nt * L * L / 2)) for nt in N_TARGETS}))
        for U in US:
            for d in DELTAS:
                jobs = [(L, nup, d, U, s) for nup in nups for s in SEEDS]
                t0 = time.time()
                with Pool(NPROC) as pool:
                    out = pool.map(run_point, jobs)
                rows = [r for chunk in out for r in chunk]
                nu = "-".join(str(x) for x in nups)
                f = f"pr_L{L}_U{U:g}_d{d:g}_nup{nu}.csv"
                with open(f, "w", newline="") as fh:
                    w = csv.writer(fh); w.writerow(COLS); w.writerows(rows)
                print(f"  L={L} U={U:g} delta={d:g}: {len(rows)} rows -> {f}"
                      f"  ({(time.time()-t0)/60:.1f} min)", flush=True)
    print("all done", flush=True)
