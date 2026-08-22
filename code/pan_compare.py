"""Pairing channels of Pan et al. (arXiv:2409.16523) on our checkerboard model.

WHY THIS EXISTS. Reference [Pan] computes the pair susceptibility and the pairing
vertex for the checkerboard Hubbard model by determinant QMC, and reports d-wave
as the leading channel. That is the closest published work to our pairing result,
so the difference has to be established rather than asserted.

WHERE THE TWO MODELS ACTUALLY COINCIDE. Their Hamiltonian has ONE diagonal per
sublattice, A along +/-(1,-1) and B along +/-(1,1), the other absent. Ours has both,
t1 -/+ delta. The two are the same Hamiltonian only where our second diagonal
switches off, delta = t1 = 0.3, verified by identical spectra to 4e-15. Their t' is
positive and ours negative in their convention, but on this lattice the
particle-hole transformation c_i -> (-1)^(x+y) c_i^dagger sends t' -> -t' and
n -> 2-n, so AT HALF FILLING the two signs are equivalent and the match is exact.
Their working density n = 0.9 maps to n = 1.1, so their doped data is not ours.

WHAT IS ADDED. Their four form factors (their Eq. 8) carry weights g_I = sqrt(3)/2
on the four NN bonds and 1/2 on the two NNN bonds:

    f_s = +1 (NN),  -1 (NNN)
    f_d = (-1)^I (NN),  0 (NNN)          <- NN d_x2-y2, no diagonal weight

Their f_d is our d_x2-y2 up to the overall sqrt(3)/2. Their f_s has no counterpart
in our basis: it mixes NN with a NEGATIVE weight on the diagonals. Their complex
d+id and d+is channels are deliberately NOT implemented. They lead only above
t'/t ~ 1.2, far from the single overlap point at 0.6, and complex form factors
would require changing the real-valued contraction in _chid_block_multi, which
risks the results already committed.

Crucially, their basis contains NO d_xy. That is the channel our polarisation
analysis says the magnetism lives in, so their calculation could not have found it.

  L=12 NUPS=72 US=3,4,4.5 DELTAS=0.3 NSEED=6 NPROC=18 python pan_compare.py
"""
import os, sys, time, csv
os.environ["OMP_NUM_THREADS"] = "1"; os.environ["MKL_NUM_THREADS"] = "1"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from multiprocessing import Pool
from cpqmc import CPMC
import checkerboard as cb

T0, T1 = -1.0, 0.3
NW   = int(os.environ.get("NW", 500))
NEQ  = int(os.environ.get("NEQ", 640))       # beta = NEQ * DT = 32, matched to Fortran
NBLK = int(os.environ.get("NBLK", 40))
BP   = int(os.environ.get("BP", 16))
DT   = float(os.environ.get("DT", 0.05))
NPROC = int(os.environ.get("NPROC", min(os.cpu_count(), 18)))
NSEED = int(os.environ.get("NSEED", 6))
L = int(os.environ.get("L", 12))
US = [float(x) for x in os.environ.get("US", "3,4,4.5").split(",")]
DELTAS = [float(x) for x in os.environ.get("DELTAS", "0.3").split(",")]
NUPS = [int(x) for x in os.environ.get("NUPS", str(L * L // 2)).split(",")]

G_NN, G_NNN = np.sqrt(3.0) / 2.0, 0.5


def pan_bond_factors(lx, ly):
    """Pan et al. Eq. (8) s-wave and d-wave form factors, with their g_I weights.

    NNN bonds are SUBLATTICE-SPECIFIC: +/-(1,-1) on A (x+y even), +/-(1,1) on B,
    which are exactly the diagonals carrying t' in their Hamiltonian, and exactly
    the surviving diagonals of our model at delta = 0.3.
    """
    n = lx * ly
    Fs = np.zeros((n, n)); Fd = np.zeros((n, n))
    def idx(x, y): return (x % lx) * ly + (y % ly)
    for x in range(lx):
        for y in range(ly):
            m = idx(x, y)
            for (j, fd) in [(idx(x + 1, y), +1.0), (idx(x - 1, y), +1.0),
                            (idx(x, y + 1), -1.0), (idx(x, y - 1), -1.0)]:
                Fs[m, j] += G_NN * 1.0
                Fd[m, j] += G_NN * fd
            diag = ((1, -1), (-1, 1)) if (x + y) % 2 == 0 else ((1, 1), (-1, -1))
            for dx, dy in diag:
                Fs[m, idx(x + dx, y + dy)] += G_NNN * (-1.0)
    return Fs, Fd


def run_point(args):
    L_, nup, delta, U, seed = args
    n = L_ * L_
    K = cb.checkerboard_hopping(L_, L_, T0, T1, -delta)
    Fs, Fd = cb.nn_bond_factors(L_, L_)
    Fdxy = cb.diag_bond_factors(L_, L_)
    Ps, Pd = pan_bond_factors(L_, L_)
    Ffac = {"son": np.eye(n), "sext": Fs, "d": Fd, "dxy": Fdxy,
            "pan_s": Ps, "pan_d": Pd}
    q = CPMC(L_, L_, nup, nup, U=U, dt=DT, nwalkers=NW, seed=seed, K=K, K_dn=None)
    r = cb.run_bp_chid_cb(q, Ffac, nequil=NEQ, nblocks=NBLK, bp=BP)
    row = dict(L=L_, nup=nup, n=2 * nup / n, delta=delta, U=U, seed=seed,
               beta_proj=NEQ * DT, tau_max=BP * DT, nw=NW)
    for c in Ffac:
        row[f"chi_{c}_full"] = float(r[f"chi_{c}"])
        row[f"chi_{c}_vertex"] = float(r[f"chi_{c}_vertex"])
    return row


if __name__ == "__main__":
    jobs = [(L, nu, d, U, s) for nu in NUPS for d in DELTAS
            for U in US for s in range(1, NSEED + 1)]
    print(f"Pan-channel comparison: L={L} deltas={DELTAS} US={US} x {NSEED} seeds "
          f"= {len(jobs)} jobs | beta={NEQ*DT:g} nw={NW}", flush=True)
    t0 = time.time(); rows = []
    with Pool(NPROC) as p:
        for i, r in enumerate(p.imap_unordered(run_point, jobs), 1):
            rows.append(r)
            lead = max(("son", "sext", "d", "dxy"), key=lambda c: r[f"chi_{c}_vertex"])
            leadp = max(("pan_s", "pan_d"), key=lambda c: r[f"chi_{c}_vertex"])
            print(f"  [{i}/{len(jobs)}] U={r['U']} d={r['delta']} s={r['seed']} "
                  f"| ours: {lead}={r[f'chi_{lead}_vertex']:.3f} "
                  f"| Pan basis: {leadp}={r[f'chi_{leadp}_vertex']:.3f} "
                  f"({(time.time()-t0)/60:.1f}m)", flush=True)
            with open("pan_compare_partial.csv", "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
                w.writeheader(); w.writerows(rows)
    out = f"pan_compare_L{L}.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {out} ({len(rows)} rows) in {(time.time()-t0)/60:.1f} min\n")
    print(f"{'U':>5}{'delta':>7} | " + "".join(f"{c:>10}" for c in
          ("son","sext","d","dxy","pan_s","pan_d")))
    for U in US:
        for d in DELTAS:
            s = [r for r in rows if r["U"] == U and r["delta"] == d]
            if not s: continue
            print(f"{U:5.1f}{d:7.1f} | " + "".join(
                f"{np.mean([r[f'chi_{c}_vertex'] for r in s]):10.4f}"
                for c in ("son","sext","d","dxy","pan_s","pan_d")))
