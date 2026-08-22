"""Momentum-domain spin polarisation Delta_tot, via a staggered pinning field.

The altermagnetic order parameter, PRB Eq. (3):

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
# Boundary phase, +1 periodic / -1 antiperiodic. checkerboard_hopping has taken
# apx/apy all along and the pairing drivers (checkerboard_apbc_halffilling.py,
# checkerboard_twist.py) already use them -- only this Delta_tot driver was
# still hard-wired to periodic, so the shell-effect test could not be run on the
# one quantity that needs it. At U=0, half filling, PERIODIC boundaries the gap
# is EXACTLY zero at every L and delta (2 to 12 states at E_F), so the
# non-interacting ground state is degenerate and which states fill is decided by
# the trial, not the physics. Antiperiodic shifts the k-mesh off the
# high-symmetry points and closes the shell -- everywhere at L=12, but NOT at
# L=10 for delta = 0.1 and 0.2, where 4 degenerate states survive both choices.
APX = int(os.environ.get("APX", 1))
APY = int(os.environ.get("APY", 1))


def stagger(L_):
    """S_ii = (-1)^(x+y) with the site ordering i = x*ly + y."""
    i = np.arange(L_ * L_)
    return (-1.0) ** ((i // L_) + (i % L_))


def _nk(G, L_, offx, offy):
    """n(k) = (1/N) sum_ij exp(-i k.(r_i - r_j)) G_ij, evaluated on the k-mesh
    k = (2pi/L)(m + off).

    off = 0 for periodic, 1/2 for antiperiodic. This offset is NOT cosmetic. The
    antiperiodic Hamiltonian is built by flipping the sign of bonds that cross
    the edge, so G is antiperiodic in the site difference and the folded-FFT
    transform (which silently assumes periodicity) evaluates the occupations on
    the unshifted mesh instead. Both meshes are complete bases so nothing blows
    up, but the antiperiodic eigenstates get smeared across the wrong grid: at
    U = 0, h = 0.2, L = 12 that inflates Delta_tot by 1.8-3.6x and, far worse,
    drives A_odd from ~0 (clean d-wave) to ~1 (reads as no d-wave at all). At
    off = 0 this reproduces the previous FFT result to machine precision, so all
    existing periodic data is unaffected."""
    N = L_ * L_
    i = np.arange(N); x, y = i // L_, i % L_
    mx = (2 * np.pi / L_) * (np.arange(L_) + offx)
    my = (2 * np.pi / L_) * (np.arange(L_) + offy)
    kx, ky = np.meshgrid(mx, my, indexing="ij")
    P = np.exp(-1j * (kx.ravel()[:, None] * x[None, :]
                      + ky.ravel()[:, None] * y[None, :]))
    nk = np.real(((P @ np.asarray(G, float)) * P.conj()).sum(1)) / N
    return nk.reshape(L_, L_)


def delta_tot(Gu, Gd, L_, apx=1, apy=1):
    """Delta_tot = sum_k |n_up(k) - n_dn(k)| on the mesh set by the boundary
    condition."""
    ox = 0.0 if apx == 1 else 0.5
    oy = 0.0 if apy == 1 else 0.5
    a = _nk(Gu, L_, ox, oy); b = _nk(Gd, L_, ox, oy)
    return float(np.abs(a - b).sum()), a, b


def run_point(args):
    L_, nup, delta, U, h, seed = args
    n = L_ * L_
    K = cb.checkerboard_hopping(L_, L_, T0, T1, -delta, apx=APX, apy=APY)
    S = np.diag(stagger(L_))
    q = CPMC(L_, L_, nup, nup, U=U, dt=DT, nwalkers=NW, seed=seed,
             K=K + h * S, K_dn=K - h * S)
    r = q.run_bp_obs(nequil=NEQ, nblocks=NBLK, bp=BP)
    Gu = np.array(r["green_up"]); Gd = np.array(r["green_dn"])
    dt_, nku, nkd = delta_tot(Gu, Gd, L_, apx=APX, apy=APY)
    # staggered real-space moment, for cross-reference
    s = stagger(L_)
    m_stag = float(np.abs((np.diag(Gu) - np.diag(Gd)) * s).sum() / n)
    row = dict(L=L_, nup=nup, n=2 * nup / n, delta=delta, U=U, h=h, seed=seed,
               apx=APX, apy=APY, delta_tot=dt_, m_stag=m_stag)
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
    # boundary phase in the filename AND in the rows: a periodic and an
    # antiperiodic run must never be mistaken for one another after the fact
    bc = "" if (APX, APY) == (1, 1) else f"_apx{APX}apy{APY}"
    out = f"polarization_L{L}{bc}.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    npz = f"dnk_L{L}{bc}.npz"
    np.savez_compressed(npz, dnk=dnk,
                        meta=np.array([[r["U"], r["delta"], r["h"], r["seed"], r["n"]]
                                       for r in rows]),
                        cols=np.array(["U", "delta", "h", "seed", "n"]))
    print(f"wrote {out} ({len(rows)} rows) and {npz} "
          f"(dn(k) grids, {dnk.shape}) in {(time.time()-t0)/60:.1f} min", flush=True)
