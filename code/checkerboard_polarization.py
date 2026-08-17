"""Momentum-domain spin polarisation Delta_tot, via a staggered pinning field.

The altermagnetic order parameter, PRL Eq. (2):

    Delta_tot = sum_k |n_up(k) - n_dn(k)|

The absolute value is INSIDE the sum, so it is nonzero even though the net
magnetisation sum_k [n_up(k) - n_dn(k)] vanishes. Compensated but spin-split:
the altermagnetic signature.

Why it needs a pinning field here. Our checkerboard anisotropy is
spin-INdependent and acts on the second-neighbour (diagonal) bonds. That is the
EMERGENT altermagnetism route of arXiv 2607.06106, where the splitting appears
only once antiferromagnetic order forms, because the ordered moments together
with the C4-related sublattices break the spin degeneracy. It is NOT the PRL's
route, where spin-dependent nearest-neighbour hopping splits the bands by hand
at U = 0.

So Delta_tot must be zero at U = 0 (correct physics), and it is also zero in our
production runs for a separate reason: they use n_up = n_dn with K_up = K_dn, so
the two spin species are identical by construction and no order can break that.

A staggered field h couples oppositely to the two sublattices,
    K_up = K + h*S,  K_dn = K - h*S,   S_ii = (-1)^(x+y),
which makes the trial state symmetry-broken (the engine builds the trial from
each spin's own matrix) and lets the emergent splitting appear. Delta_tot is
then extrapolated to h -> 0. A nonzero intercept is spontaneous altermagnetism;
an intercept consistent with zero means the splitting was entirely field-induced.

CAUTION: the constrained path is defined by the trial, so a symmetry-broken
trial biases the walkers toward the broken state. The h -> 0 extrapolation is
mandatory, not optional, and the result should be read against the finite-size
scaling, which currently shows S(q*) flat in N (no long-range order).

  L=12 US=4 DELTAS=0.4 HS=0.5,0.3,0.2,0.1,0.05 NUPS=72 NSEED=6 NPROC=30 \
      python checkerboard_polarization.py
"""
import os, sys, time, csv
os.environ["OMP_NUM_THREADS"] = "1"; os.environ["MKL_NUM_THREADS"] = "1"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from multiprocessing import Pool
from cpqmc import CPMC
import checkerboard as cb

T0, T1 = -1.0, 0.3
NW, NEQ, NBLK, BP, DT = 160, 60, 30, 12, 0.05
NPROC = int(os.environ.get("NPROC", min(os.cpu_count(), 30)))
NSEED = int(os.environ.get("NSEED", 6))
L = int(os.environ.get("L", 12))
US = [float(x) for x in os.environ.get("US", "2,4,6,8").split(",")]
DELTAS = [float(x) for x in os.environ.get("DELTAS", "0,0.2,0.4,0.6").split(",")]
HS = [float(x) for x in os.environ.get("HS", "0.5,0.3,0.2,0.1,0.05").split(",")]
NUPS = [int(x) for x in os.environ.get("NUPS", str(L * L // 2)).split(",")]


def stagger(L_):
    """S_ii = (-1)^(x+y) with the site ordering i = x*ly + y."""
    i = np.arange(L_ * L_)
    return (-1.0) ** ((i // L_) + (i % L_))


def delta_tot(Gu, Gd, L_):
    """Delta_tot = sum_k |n_up(k) - n_dn(k)| from the real-space Green functions.
    n_s(k) = (1/N) sum_ij exp(-i k.(r_i - r_j)) G^s_ij, obtained by a 2D FFT of
    G reshaped onto the lattice."""
    N = L_ * L_
    idx = np.arange(N)
    x, y = idx // L_, idx % L_
    out = []
    for G in (Gu, Gd):
        # accumulate G_ij into the translation-averaged correlator g(dx,dy)
        g = np.zeros((L_, L_))
        dx = (x[:, None] - x[None, :]) % L_
        dy = (y[:, None] - y[None, :]) % L_
        np.add.at(g, (dx.ravel(), dy.ravel()), np.asarray(G).ravel())
        g /= N
        out.append(np.real(np.fft.fft2(g)))
    return float(np.abs(out[0] - out[1]).sum()), out[0], out[1]


def run_point(args):
    L_, nup, delta, U, h, seed = args
    n = L_ * L_
    K = cb.checkerboard_hopping(L_, L_, T0, T1, -delta)
    S = np.diag(stagger(L_))
    q = CPMC(L_, L_, nup, nup, U=U, dt=DT, nwalkers=NW, seed=seed,
             K=K + h * S, K_dn=K - h * S)
    r = q.run_bp_obs(nequil=NEQ, nblocks=NBLK, bp=BP)
    Gu = np.array(r["green_up"]); Gd = np.array(r["green_dn"])
    dt_, nku, nkd = delta_tot(Gu, Gd, L_)
    # staggered real-space moment, for cross-reference
    s = stagger(L_)
    m_stag = float(np.abs((np.diag(Gu) - np.diag(Gd)) * s).sum() / n)
    row = dict(L=L_, nup=nup, n=2 * nup / n, delta=delta, U=U, h=h, seed=seed,
               delta_tot=dt_, m_stag=m_stag)
    # KEEP the momentum-resolved occupations. delta_tot() computes n_up(k) and
    # n_dn(k) anyway, and without them only the scalar Delta_tot survives -- the
    # Delta n(k) map over the Brillouin zone, which is what the PRL Fig. 1 and
    # our own schematic actually show, would be lost.
    return row, (nku - nkd)


if __name__ == "__main__":
    jobs = [(L, nu, d, U, h, s) for nu in NUPS for d in DELTAS
            for U in US for h in HS for s in range(1, NSEED + 1)]
    print(f"Delta_tot pinning-field scan: L={L} nups={NUPS} US={US} "
          f"deltas={DELTAS} h={HS} x {NSEED} seeds = {len(jobs)} jobs", flush=True)
    t0 = time.time()
    with Pool(NPROC) as p:
        res = p.map(run_point, jobs)             # list of (row, dn(k))
    rows = [r for r, _ in res]
    dnk = np.array([g for _, g in res])          # (njobs, L, L)
    out = f"polarization_L{L}.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    npz = f"dnk_L{L}.npz"
    np.savez_compressed(npz, dnk=dnk,
                        meta=np.array([[r["U"], r["delta"], r["h"], r["seed"], r["n"]]
                                       for r in rows]),
                        cols=np.array(["U", "delta", "h", "seed", "n"]))
    print(f"wrote {out} ({len(rows)} rows) and {npz} "
          f"(dn(k) grids, {dnk.shape}) in {(time.time()-t0)/60:.1f} min", flush=True)
