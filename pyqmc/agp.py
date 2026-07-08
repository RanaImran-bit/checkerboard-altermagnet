#!/usr/bin/env python3
"""AGP / number-projected-BCS geminal trial wavefunction for CP-AFQMC.

The d-wave pairing constrained-path bias (ED +23.18 vs CP-AFQMC ~ -2 on the
4x2/t1=0.3/U=4 gate) is the project's key open problem: a free-electron trial has
the wrong NODES for the 4-point pairing channel, and constraint release fixes only
the ket, not the bra (see docs/single_band_altermagnet_report.md). The fix is a
pairing-aware trial. This module supplies the rigorous overlap of a paired trial
with a number-conserving Slater walker -- the piece the sister Nambu/BdG effort
(thkaitools-cell/pph-qmc, report sec.6) was missing.

THE STATE (antisymmetrized geminal power = number-projected BCS):
    |Psi_AGP> = ( sum_{ij} F_ij c^+_{i up} c^+_{j dn} )^N |0> ,   N = nup = ndn,
F (n x n) the geminal matrix. For a Slater walker |Phi> with up/dn orbital
matrices A (n x N), B (n x N):

    <Psi_AGP | Phi> = det(A^T F B)          (real-valued / theta=0 convention)

The N-column determinant automatically projects |BCS> onto the N-pair sector.

MIXED CONTRACTIONS (the subtle part -- validated to 1e-15 vs ED in validate_agp.py):
the projected-AGP is NOT Gaussian, so a normal-only Wick is WRONG. The fix is the
identity that, because H conserves particle number and the walker ket |Phi> has
fixed (N,N), <Psi_AGP|H|Phi>/<Psi_AGP|Phi> = <BCS|H|Phi>/<BCS|Phi> with the
UNPROJECTED Gaussian |BCS> = exp(sum F c^+ c^+)|0>. Generalized Wick on that
Gaussian gives NORMAL densities
    M = A^T F B
    Gu[i,j] = <c^+_{i up} c_{j up}> = (F B M^{-1} A^T)_{ij}
    Gd[i,j] = <c^+_{i dn} c_{j dn}> = (F^T A M^{-T} B^T)_{ij}
PLUS the nonzero ANOMALOUS (pairing) tensors
    Kd[i,j] = <c^+_{i up} c^+_{j dn}> = (F (I - Gd^T))_{ij}
    Ka[i,j] = <c_{i up}   c_{j dn}>   = -(A M^{-T} B^T)_{ij}
The cross-spin density-density term then carries an anomalous correction:
    <n_{a up} n_{b dn}> = Gu[a,a] Gd[b,b] - Kd[a,b] Ka[a,b]
(same-spin terms keep the ordinary Wick exchange; singlet pairing only couples
up-dn). This -Kd*Ka piece is exactly what a normal-only estimator (and the sister
pph-qmc Nambu prototype, report sec.6) was missing. (Built-in sanity: tr(Gu)=tr(Gd)=N.)

GEMINAL FROM A d-WAVE BdG MEAN FIELD (cf. pph-qmc hf_bdg_twist_prototype.py):
Nambu spinor Psi=(c_up, c^+_dn),  H_BdG = [[h_up - mu, Delta], [Delta^dag, -h_dn^* + mu]].
The positive-energy eigenvectors stack as columns (U; V) (each n x n); the
quasiparticle-vacuum condition (U^dag c_up + V^dag c^+_dn)|GS> = 0 gives
    F = -(U^dag)^{-1} V^dag .
Delta is the d-wave bond field (+x bonds, -y bonds) with amplitude lambda_d.

This module is intentionally standalone (pure functions on numpy arrays); the
validated cpqmc.py core is untouched so its bit-for-bit regression gate is unaffected.
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from cpqmc import am_hopping            # spin-dependent anisotropic + NNN hopping


# ---------------------------------------------------------------- d-wave BdG
def dwave_delta(lx, ly, lam):
    """Real d-wave singlet pairing field Delta_ij = lam * f_d(bond), f_d = +1 on
    x bonds, -1 on y bonds (symmetric in i<->j, as a singlet pair field must be).
    Site index i = x*ly + y to match am_hopping / cpqmc."""
    n = lx * ly
    D = np.zeros((n, n))
    def idx(x, y): return (x % lx) * ly + (y % ly)
    for x in range(lx):
        for y in range(ly):
            i = idx(x, y)
            for (dx, dy, fd) in ((1, 0, +1.0), (-1, 0, +1.0), (0, 1, -1.0), (0, -1, -1.0)):
                j = idx(x + dx, y + dy)
                D[i, j] += 0.5 * lam * fd     # 0.5 b/c each bond hit from both ends
    return D


def bdg_geminal(h_up, h_dn, delta, mu=0.0, reg=1e-10):
    """Assemble the BdG matrix for Nambu Psi=(c_up, c^+_dn), diagonalize, and
    return the geminal F = -(U^dag)^{-1} V^dag from the positive-energy block.

      H_BdG = [[h_up - mu,      Delta        ],
               [Delta^dag,  -h_dn^* + mu      ]]

    h_up, h_dn are the one-body (hopping) matrices; delta the pairing field.
    Returns (F, info) with info carrying the spectrum + diagnostics."""
    n = h_up.shape[0]
    I = np.eye(n)
    H = np.block([[h_up - mu * I, delta],
                  [delta.conj().T, -h_dn.conj() + mu * I]])
    herm = float(np.max(np.abs(H - H.conj().T)))
    w, V = np.linalg.eigh(H)
    # positive-energy eigenvectors -> quasiparticle creation; their (U;V) blocks
    # define the quasiparticle-vacuum geminal. Take the n largest eigenvalues.
    order = np.argsort(w)[::-1]            # descending
    pos = order[:n]
    Wp = V[:, pos]
    U = Wp[:n, :]
    Vb = Wp[n:, :]
    # F = -(U^dag)^{-1} V^dag.  For a NODAL (d-wave) gap the U block is singular
    # (node/deep orbitals are unpaired, v/u -> inf), so use a Tikhonov-regularized
    # pseudo-inverse: those directions get a large-but-finite (near-frozen-core)
    # pairing amplitude. The trial only needs the right NODES, not the exact BCS
    # ground state, so this regularization is physically harmless.
    Udag = U.conj().T
    UU = Udag @ U + reg * np.eye(n)
    F = -np.linalg.solve(UU, Udag @ Vb.conj().T)
    info = dict(evals=w, bdg_herm_err=herm, U=U, V=Vb, reg=reg,
                min_abs_eval=float(np.min(np.abs(w))),
                Ucond=float(np.linalg.cond(U)),
                Fnorm=float(np.linalg.norm(F)))
    return F, info


def dwave_geminal(lx, ly, nup, ndn, t0=1.0, tam=0.0, t1=0.0, lam=0.1, mu=0.0,
                  h_up=None, h_dn=None):
    """Convenience: build the bare anisotropic+NNN band (am_hopping), add a d-wave
    pairing field of amplitude `lam`, and return the geminal F. Pass h_up/h_dn to
    override the one-body part with a (UHF) mean field instead of the bare band."""
    if h_up is None or h_dn is None:
        Ku, Kd = am_hopping(lx, ly, t0, tam, t1)
        h_up = Ku if h_up is None else h_up
        h_dn = Kd if h_dn is None else h_dn
    delta = dwave_delta(lx, ly, lam)
    F, info = bdg_geminal(h_up, h_dn, delta, mu=mu)
    return F, info


def augmented_geminal(phi_up, phi_dn, Fd, eta):
    """A WELL-CONDITIONED d-wave-correlated geminal built as the free-electron
    determinant PLUS a bounded d-wave pairing admixture:

        F = phi_up phi_dn^T + eta * Fd

    At eta=0 this reproduces the free trial EXACTLY: det(A^T F B) =
    det(A^T phi_up) det(phi_dn^T B) = <phi_up|A><phi_dn|B> (factorized Slater
    overlap), so CP-AFQMC recovers the validated fixed-trial result. Increasing eta
    smoothly turns on d-wave pairing correlations without the nodal v/u -> inf
    singularity of the bare BdG geminal. phi_up/phi_dn are n x N orbital matrices;
    Fd the d-wave bond form-factor matrix."""
    return phi_up @ phi_dn.T + eta * Fd


# ---------------------------------------------------------------- AGP overlap / GF
def agp_overlap(F, A, B):
    """<Psi_AGP|Phi> = det(A^T F B)  for walker up/dn orbital matrices A,B (n x N).
    Real-valued (theta=0) convention; returns a scalar (may be negative)."""
    M = A.T @ F @ B
    sign, logdet = np.linalg.slogdet(M)
    return sign * np.exp(logdet)


def agp_green(F, A, B):
    """Mixed one-body density matrices for the AGP bra / Slater ket, fixed (N,N):
        Gu[i,j] = <c^+_{i up} c_{j up}> = (F B M^{-1} A^T)_{ij}
        Gd[i,j] = <c^+_{i dn} c_{j dn}> = (F^T A M^{-T} B^T)_{ij}
    M = A^T F B. Returns (Gu, Gd) as n x n arrays in the SAME index convention the
    cpqmc Estimators use after their `.T` (Gu[i,j] = <c^+_i c_j>)."""
    M = A.T @ F @ B
    Minv = np.linalg.inv(M)
    Gu = F @ B @ Minv @ A.T
    Gd = (F.T @ A @ Minv.T @ B.T)
    return Gu, Gd


def agp_overlap_green(F, A, B):
    """Overlap + both Green's functions sharing one M factorization."""
    M = A.T @ F @ B
    sign, logdet = np.linalg.slogdet(M)
    O = sign * np.exp(logdet)
    Minv = np.linalg.inv(M)
    Gu = F @ B @ Minv @ A.T
    Gd = F.T @ A @ Minv.T @ B.T
    return O, Gu, Gd


def agp_contractions(F, A, B):
    """All mixed contractions of the AGP bra / Slater ket needed for the energy and
    pairing observables. Returns dict with overlap O and matrices (n x n):
        Gu[i,j] = <c^+_{i up} c_{j up}>          Gd[i,j] = <c^+_{i dn} c_{j dn}>
        Kd[i,j] = <c^+_{i up} c^+_{j dn}>        Ka[i,j] = <c_{i up} c_{j dn}>
    """
    M = A.T @ F @ B
    sign, logdet = np.linalg.slogdet(M)
    O = sign * np.exp(logdet)
    Minv = np.linalg.inv(M)
    Gu = F @ B @ Minv @ A.T
    Gd = F.T @ A @ Minv.T @ B.T
    Kd = F @ (np.eye(F.shape[0]) - Gd.T)         # <c^+_up c^+_dn>
    Ka = -A @ Minv.T @ B.T                       # <c_up c_dn>
    return dict(O=O, Gu=Gu, Gd=Gd, Kd=Kd, Ka=Ka)


def agp_local_energy(model, F, A, B):
    """Rigorous mixed local energy E_L = <Psi_AGP|H|Phi>/<Psi_AGP|Phi> for the AGP
    trial, using model.K / model.K_dn (one-body) and model.pot_terms
    (V n_{a,sa} n_{b,sb}). Cross-spin density-density terms carry the anomalous
    correction -Kd[a,b]*Ka[a,b]; same-spin keep the ordinary Wick exchange. This is
    the generalized-Slater/HFB local energy (validated vs ED in validate_agp.py)."""
    c = agp_contractions(F, A, B)
    Gu, Gd, Kd, Ka = c["Gu"], c["Gd"], c["Kd"], c["Ka"]
    e = float(np.sum(model.K * Gu) + np.sum(model.K_dn * Gd))   # Gu[i,j]=<c^+_i c_j>
    for (a, sa, b, sb, V) in model.pot_terms:
        if sa != sb:                              # up-dn: factorized minus anomalous
            e += V * (Gu[a, a] * Gd[b, b] - Kd[a, b] * Ka[a, b]) if (sa, sb) == (0, 1) \
                 else V * (Gd[a, a] * Gu[b, b] - Kd[b, a] * Ka[b, a])
        else:                                     # same spin: ordinary Wick exchange
            Gs = Gu if sa == 0 else Gd
            e += V * (Gs[a, a] * Gs[b, b] - Gs[a, b] * Gs[b, a])
    return e
