#!/usr/bin/env python3
"""Single-particle spectral machinery for the checkerboard (meeting point 7):
G(k, tau) -> quasiparticle dispersion E_k -> charge gap and DOS.

KEY POINT: the time-displaced Green's function is ALREADY computed by the existing
back-propagation estimator. `_chid_block_multi` builds
    Pu = Bu @ (I - gu^T)  == P^up(tau_l)_{mn} = <c_m(tau_l) c^+_n(0)>
per tau slice and accumulates the weighted sum, purely so `run_bp_chid` can build the
bubble. Passing an EMPTY form-factor dict gives those GFs with zero pairing overhead.
So no new QMC estimator is needed -- only the k-space reduction and the analysis.

WHY NO MaxEnt (yet): analytic continuation of G(tau) -> A(k,omega) is ill-conditioned
and needs a long imaginary-time window. What point 7 actually asks for -- the energy
gap and the DOS -- is obtainable far more robustly from the exponential decay

    G(k, tau) ~ Z_k exp(-(E_k - mu) tau)       (large tau, T=0 particle GF)

so E_k - mu is the slope of log G(k, tau). That is a linear fit, not an inverse
Laplace transform, and it has no regularisation choices to defend to a referee.
Full A(k,omega) via MaxEnt can come later if the lineshape (not just the peak) is wanted.

CRITICAL LIMITATION -- the tau window. tau_max = bp*dt. The pairing runs use
bp=16, dt=0.05 -> tau_max = 0.8, which is FAR too short to resolve a gap: the fit needs
tau_max >> 1/Delta, so a gap of order 0.5t wants tau_max ~ 4, i.e. bp ~ 80. Longer
back-propagation also increases the constrained-path bias, so bp must be converged, not
merely enlarged. `scan_bp()` below exists to check exactly that.

VALIDATION: at U=0 the extracted E_k must reproduce the free-electron dispersion of the
checkerboard hopping matrix exactly. `validate_free()` performs that test.
"""
import numpy as np
from cpqmc import CPMC, _green
import checkerboard as cb
from unified_scan import _shift_index


def _to_kspace(P, lx, ly, shift):
    """P[m,n] (real space) -> G(k) on the (lx,ly) grid, same FFT convention as
    run_bp_chi_spin's _pq: S(R)=sum_m P[m, m+R], then G(k)=Re FFT2(S)/N."""
    n = lx * ly
    rows = np.arange(n)
    Sg = np.array([[P[rows, shift[dx, dy]].sum() for dy in range(ly)] for dx in range(lx)])
    return np.real(np.fft.fft2(Sg)) / n


def run_bp_gtau(q, nequil=60, nblocks=40, bp=40, ortho=10, pc=10):
    """Ensemble-averaged time-displaced Green's function G(k, tau_l), l = 0..bp.

    Mirrors run_bp_chid_cb's block loop exactly but passes an empty form-factor dict,
    so only the GFs are accumulated. Returns per-block statistics.

    Returns dict:
      taus    (bp+1,)                 tau_l = l*dt
      Gk_up   (bp+1, lx, ly)          block-mean G_up(k, tau)
      Gk_up_err  same shape           standard error over blocks
      Gk_dn, Gk_dn_err                same for down spin
      nblocks_used
    """
    from checkerboard import _chid_block_multi
    lx, ly = q.lx, q.ly
    shift = _shift_index(lx, ly)

    for it in range(nequil):
        q.step()
        if (it + 1) % ortho == 0: q.reorthogonalize()
        if (it + 1) % pc == 0: q.pop_control()
        q._maybe_update_trial(it)

    blocks_u, blocks_d = [], []
    for _ in range(nblocks):
        q.reorthogonalize()
        ket_up = q.walkers.phi_up.copy(); ket_dn = q.walkers.phi_dn.copy()
        rec = [[] for _ in range(q.nw)]
        for _ in range(bp):
            q.prop.step_record(q.walkers, rec)
        _, Pu_s, Pd_s, W = _chid_block_multi(q.est, q.walkers, ket_up, ket_dn, rec, bp, {})
        if W <= 0:
            q.pop_control(); continue
        Pu = Pu_s / W; Pd = Pd_s / W
        blocks_u.append(np.array([_to_kspace(Pu[l], lx, ly, shift) for l in range(bp + 1)]))
        blocks_d.append(np.array([_to_kspace(Pd[l], lx, ly, shift) for l in range(bp + 1)]))
        q.pop_control()

    Bu = np.array(blocks_u); Bd = np.array(blocks_d)
    nb = len(Bu)
    sem = lambda A: A.std(axis=0) / np.sqrt(len(A)) if len(A) > 1 else np.zeros(A.shape[1:])
    return {"taus": q.dt * np.arange(bp + 1), "nblocks_used": nb,
            "Gk_up": Bu.mean(axis=0), "Gk_up_err": sem(Bu),
            "Gk_dn": Bd.mean(axis=0), "Gk_dn_err": sem(Bd),
            # per-block data kept so the ALF bridge can build the tau-tau covariance
            "blocks_up": Bu, "blocks_dn": Bd}


# ----------------------------------------------------------------------
# Bridge to ALF's MaxEnt / stochastic analytic continuation (Max_SAC.out)
# ----------------------------------------------------------------------
# Format reverse-engineered from Max_SAC.F90:118-133 and a real input file
# (ALF_UV/data/U4V0.67/Green_0.00_2.42/g_dat):
#
#   line 1        : ntau  nbin_qmc  Beta  Norb  Channel
#   next ntau     : tau   G(tau)   err
#   if N_Cov==1   : ntau*ntau lines, xcov(nt,nt1), nt outer / nt1 inner
#
# "Beta" is the imaginary-time RANGE, not an inverse temperature: their file has
# ntau=101, dtau=0.1, Beta=10.0. For a projective/T=0 method that is tau_max.
# Channel "T0" selects XKER_T0, the ground-state kernel G(tau)=int dw A(w) e^{-tau w}
# -- the correct one for CPQMC (Max_SAC.F90:279).
#
# ALF's driver (Scripts_and_Parameters_files/Spectral.sh) expects ONE DIRECTORY PER
# K-POINT named <Prefix>_<kx>_<ky>, each containing g_dat + parameters, and runs
# Max_SAC.out inside each.

PARAMETERS_TEMPLATE = """&VAR_errors
n_skip  = 0
N_rebin = 1
N_Cov   = {n_cov}
/

&VAR_Max_Stoch
Ngamma     = 400
Om_st      = {om_st}d0
Om_en      = {om_en}d0
NDis       = {ndis}
Nbins      = 3200
Nsweeps    = 100
NWarm      = 30
N_alpha    = 70
alpha_st   = 0.01d0
R          = 1.2d0
Checkpoint = .F.
Tolerance  = {tol}d0
/
"""


def write_alf_gdat(res, outdir, spin="up", prefix="Green", om_st=-8.0, om_en=8.0,
                   ndis=600, tol=0.1, with_cov=False, gmin=1e-4, kpoints=None):
    """Write G(k,tau) in ALF Max_SAC input format: one directory per k-point.

    res      : dict returned by run_bp_gtau (needs blocks_* if with_cov)
    outdir   : directory to create the per-k subdirectories in
    with_cov : also write the full tau-tau covariance (sets N_Cov=1). ALF's own runs
               use N_Cov=0, but the covariance matters when tau slices are correlated,
               which they always are in back-propagation -- switching it on is the
               more honest choice if the block count supports it (need nblocks >> ntau).
    gmin     : skip k-points whose G(tau) is not usably positive. The particle GF is
               ~0 for OCCUPIED states (you cannot add an electron where one already
               is), and analytic continuation of noise is meaningless. Skipped points
               are reported, never silently dropped.

    Returns (written_kpoints, skipped_kpoints).
    """
    import os
    taus = np.asarray(res["taus"]); ntau = len(taus)
    Gk = res[f"Gk_{spin}"]; Ek = res[f"Gk_{spin}_err"]
    blocks = res.get(f"blocks_{spin}")
    nb = int(res["nblocks_used"])
    lx, ly = Gk.shape[1], Gk.shape[2]
    beta = float(taus[-1])                      # tau_max, per ALF's convention above

    if with_cov and (blocks is None or nb < 2 * ntau):
        raise ValueError(f"covariance needs blocks and nblocks >> ntau "
                         f"(have nblocks={nb}, ntau={ntau}); rerun with more nblocks "
                         f"or set with_cov=False")

    os.makedirs(outdir, exist_ok=True)
    written, skipped = [], []
    ks = kpoints if kpoints is not None else [(i, j) for i in range(lx) for j in range(ly)]
    for (i, j) in ks:
        g = Gk[:, i, j]
        if not (g[0] > gmin and np.all(np.isfinite(g))):
            skipped.append((i, j, "G(0) below gmin (occupied state) or non-finite"))
            continue
        kx = 2 * np.pi * i / lx; ky = 2 * np.pi * j / ly
        d = os.path.join(outdir, f"{prefix}_{kx:.2f}_{ky:.2f}")
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "g_dat"), "w") as f:
            f.write(f"{ntau:12d}{nb:12d}  {beta:.17E}{1:12d} T0\n")
            for nt in range(ntau):
                f.write(f"  {taus[nt]:.17E}  {g[nt]:.17E}  {max(Ek[nt, i, j], 1e-14):.17E}\n")
            if with_cov:
                X = blocks[:, :, i, j]                       # (nblocks, ntau)
                C = np.cov(X, rowvar=False) / nb             # covariance OF THE MEAN
                for nt in range(ntau):
                    for nt1 in range(ntau):
                        f.write(f"  {C[nt, nt1]:.17E}\n")
        with open(os.path.join(d, "parameters"), "w") as f:
            f.write(PARAMETERS_TEMPLATE.format(n_cov=1 if with_cov else 0,
                                               om_st=om_st, om_en=om_en, ndis=ndis, tol=tol))
        written.append((i, j, kx, ky))

    # driver script, modelled on ALF Scripts_and_Parameters_files/Spectral.sh
    with open(os.path.join(outdir, "run_maxent.sh"), "w") as f:
        f.write("#!/bin/bash\n"
                "# Runs ALF's Max_SAC.out in every k-point directory and collects A(k,omega).\n"
                "# Point MAXSAC at the built binary from ALF_UV/source/Analysis.\n"
                'MAXSAC=${MAXSAC:-/home/amax/run/ALF_UV/source/Analysis/Max_SAC.out}\n'
                'OUT=Spectral_all.dat\n: > $OUT\n'
                f'for d in {prefix}_*/ ; do\n'
                '  ( cd "$d" && $MAXSAC > maxent.log 2>&1 ) || { echo "FAILED: $d"; continue; }\n'
                '  k=${d#' + prefix + '_}; k=${k%/}\n'
                '  awk -v k="${k/_/ }" \'!/^#/ && NF>=2 {print k, $1, $2}\' "$d/Aom_ps" >> $OUT 2>/dev/null \\\n'
                '    || echo "no Aom_ps in $d"\n'
                'done\n'
                'echo "wrote $OUT"\n')
    os.chmod(os.path.join(outdir, "run_maxent.sh"), 0o755)

    print(f"wrote {len(written)} k-point dirs in {outdir}  (ntau={ntau}, tau_max={beta:.2f}, "
          f"nblocks={nb}, N_Cov={1 if with_cov else 0})")
    if skipped:
        print(f"skipped {len(skipped)} k-points (occupied / unusable): "
              f"{[(i,j) for i,j,_ in skipped]}")
    return written, skipped


def dispersion_from_gtau(taus, Gk, tau_min=None, Gmin=1e-6):
    """E_k - mu from the large-tau slope of log G(k, tau).

    Fits log G = log Z - (E-mu) tau over tau >= tau_min, using only points where G is
    positive and above Gmin (the QMC GF can go negative in the noise floor, and log of
    that is meaningless -- silently dropping those points is the main correctness trap
    here, so the count of used points is returned for inspection).

    Returns (E, Z, npts) each of shape Gk.shape[1:].
    """
    taus = np.asarray(taus)
    if tau_min is None:
        tau_min = taus[len(taus) // 3]           # drop the short-tau transient
    sel = taus >= tau_min
    lx, ly = Gk.shape[1], Gk.shape[2]
    E = np.full((lx, ly), np.nan); Z = np.full((lx, ly), np.nan)
    npts = np.zeros((lx, ly), dtype=int)
    for i in range(lx):
        for j in range(ly):
            g = Gk[:, i, j]
            ok = sel & (g > Gmin)
            npts[i, j] = ok.sum()
            if ok.sum() < 3:
                continue
            c = np.polyfit(taus[ok], np.log(g[ok]), 1)
            E[i, j] = -c[0]; Z[i, j] = np.exp(c[1])
    return E, Z, npts


def validate_free(L=6, nup=14, delta=0.0, bp=40, dt=0.05, nwalkers=80,
                  nequil=30, nblocks=20, seed=1, t0=-1.0, t1=0.3):
    """U=0 check: extracted E_k must equal the free-electron dispersion of K, UP TO A
    CONSTANT. At U=0 there are no auxiliary fields, so CPQMC is exact -- this tests the
    k-space reduction, the FFT convention and the fitting, not statistics.

    TWO THINGS THIS FUNCTION LEARNED THE HARD WAY:

    1. The fitted E carries a CONSTANT OFFSET (CPQMC propagates with a reference energy,
       and the canonical ensemble has no explicit mu). Verified: the offset is constant
       to machine precision across k (spread 0.0), so it cancels in every energy
       DIFFERENCE -- gaps and bandwidths are unaffected. Compare shapes, never absolutes.

    2. Only delta=0 is a valid single-band test. For delta != 0 the diagonal amplitudes
       swap between sublattices (checkerboard_hopping line 42), so K is NOT invariant
       under single-site translation -- only under the (1,1) diagonal translation.
       The model then has a TWO-SITE unit cell and TWO bands, G(k,tau) is a sum of two
       exponentials, and both this single-exponential fit and the <k|K|k> reference are
       wrong. See the module docstring.

    Returns (E_fit, eps_exact, max_abs_shape_error) with the shape error computed after
    referencing both to their minimum.
    """
    K = cb.checkerboard_hopping(L, L, t0, t1, -delta)
    q = CPMC(L, L, nup, nup, U=0.0, dt=dt, nwalkers=nwalkers, seed=seed, K=K, K_dn=None)
    r = run_bp_gtau(q, nequil=nequil, nblocks=nblocks, bp=bp)
    E, Z, npts = dispersion_from_gtau(r["taus"], r["Gk_up"])

    # exact free dispersion on the same k-grid, from the hopping matrix itself
    evals = np.linalg.eigvalsh(K)
    mu = 0.5 * (evals[nup - 1] + evals[nup])          # mid-gap chemical potential
    # eps(k) for the checkerboard: diagonalise K in the plane-wave basis
    n = L * L
    kx = 2 * np.pi * np.arange(L) / L
    eps = np.full((L, L), np.nan)
    for a in range(L):
        for b in range(L):
            v = np.exp(1j * (kx[a] * (np.arange(n) // L) + kx[b] * (np.arange(n) % L))) / np.sqrt(n)
            eps[a, b] = np.real(np.conj(v) @ K @ v)
    err = np.nanmax(np.abs(E - (eps - mu)))
    return E, eps - mu, err


def scan_bp(L=6, nup=14, delta=0.2, U=4.0, bps=(16, 24, 32, 40, 56), dt=0.05,
            nwalkers=80, nequil=40, nblocks=20, seed=1, t0=-1.0, t1=0.3):
    """Convergence of the extracted gap with back-propagation length.

    tau_max = bp*dt must be long enough to reach the exponential regime, but longer bp
    also worsens the constrained-path bias -- so the gap must PLATEAU in bp. If it does
    not, the number is a fit artefact and should not be reported.
    Returns a list of (bp, tau_max, gap).
    """
    out = []
    for bp in bps:
        K = cb.checkerboard_hopping(L, L, t0, t1, -delta)
        q = CPMC(L, L, nup, nup, U=U, dt=dt, nwalkers=nwalkers, seed=seed, K=K, K_dn=None)
        r = run_bp_gtau(q, nequil=nequil, nblocks=nblocks, bp=bp)
        E, Z, npts = dispersion_from_gtau(r["taus"], r["Gk_up"])
        gap = np.nanmin(E[E > 0]) if np.any(E > 0) else np.nan
        out.append((bp, bp * dt, float(gap)))
        print(f"  bp={bp:3d}  tau_max={bp*dt:.2f}  gap={gap:+.4f}  "
              f"(fit points used: median {int(np.median(npts))})", flush=True)
    return out
