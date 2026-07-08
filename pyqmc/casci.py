#!/usr/bin/env python3
"""Phase 6 step (b): scalable multi-determinant trial via self-consistent
natural-orbital CASCI (a mini-CASSCF, NO ED).

This module provides the deterministic CI engine that step (a)
(validate_casci_ed.py) proved the CPMC needs: a GOOD multi-determinant trial
removes the constrained-path bias, and walkers are a poor basis, so we build the
expansion by diagonalizing a small complete-active-space (CAS) Hamiltonian in a
rotated orbital basis.

Unlike validate_casci_ed.ci_ground_state (which assumes the interaction is
DIAGONAL density-density in the SITE basis), this works in an ARBITRARY
orthonormal orbital basis W: the Hubbard interaction becomes a full rank-4 tensor
    U_pqrs = U sum_i W_ip W_iq W_ir W_is        (up-down structure)
and the one-body part h_pq = (W^T K W)_pq is no longer diagonal. We evaluate
<D'|H|D> by Slater-Condon for a general one-body + the (up x down) two-body
operator, in OUR sorted-column determinant convention (same JW sign as
validate_casci_ed), so the resulting determinants map cleanly to site-basis
Slater matrices W[:, occ] for TrialWF.set_multidet.

Self-validating gate: CASCI is basis-invariant, so for ANY orthonormal W the
full-space (no frozen core) E0 must equal the exact ED ground-state energy.

    source tools/env.sh
    python pyqmc/casci.py --lx 2 --ly 2 --nup 2 --ndn 2 --U 4
"""
from __future__ import annotations
import argparse, os, sys
from itertools import combinations
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "ed"))
from validate_casci_ed import _hop_sign, det_from_config   # noqa: E402


def transform_hubbard(K, U, W):
    """Rotate the single-band Hubbard Hamiltonian into the orthonormal orbital
    basis whose columns are W (n_site x n_orb): one-body h = W^T K W, two-body
    U_pqrs = U sum_i W_ip W_iq W_ir W_is (couples up density a^+_p a_q to down
    a^+_r a_s)."""
    W = np.asarray(W, float)
    h1 = W.T @ K @ W
    Upqrs = U * np.einsum("ip,iq,ir,is->pqrs", W, W, W, W, optimize=True) if U != 0 \
        else np.zeros((W.shape[1],) * 4)
    return h1, Upqrs


def _spin_trans(occ_bra, occ_ket):
    """Classify the single-spin one-body transition <occ_bra| a^+_p a_q |occ_ket>.
    Returns ('diag', occ) if the occupations are equal (then a^+_p a_q is nonzero
    only for p=q in occ), ('single', p, q, sign) if they differ by exactly one
    orbital (ket has q, bra has p), or None otherwise."""
    sb, sk = set(occ_bra), set(occ_ket)
    if sb == sk:
        return ("diag", occ_ket)
    only_bra = sb - sk
    only_ket = sk - sb
    if len(only_bra) == 1 and len(only_ket) == 1:
        p = next(iter(only_bra)); q = next(iter(only_ket))
        return ("single", p, q, _hop_sign(occ_ket, p, q))
    return None


def _one_body_me(trans, h1):
    """One-body matrix element for one spin given the transition classification."""
    if trans is None:
        return 0.0
    if trans[0] == "diag":
        return sum(h1[b, b] for b in trans[1])
    _, p, q, sgn = trans
    return h1[p, q] * sgn


def _two_body_me(tu, td, Upqrs):
    """<D'|H2|D> for H2 = sum_pqrs U_pqrs a^+_p,up a_q,up a^+_r,dn a_s,dn, given the
    up/down transition classifications tu, td. Factorizes into an up one-body
    transition density times a down one-body transition density."""
    if tu is None or td is None:
        return 0.0
    # up factor: list of (p, q, weight) with <bra_up|a^+_p a_q|ket_up>
    if tu[0] == "diag":
        up = [(b, b, 1.0) for b in tu[1]]
    else:
        _, p, q, s = tu; up = [(p, q, s)]
    if td[0] == "diag":
        dn = [(b, b, 1.0) for b in td[1]]
    else:
        _, r, ss, s2 = td; dn = [(r, ss, s2)]
    val = 0.0
    for (p, q, su) in up:
        for (r, s, sd) in dn:
            val += Upqrs[p, q, r, s] * su * sd
    return val


def casci_solve(h1, Upqrs, n, nup, ndn, active=None, core=None, ecore=0.0):
    """Diagonalize the CAS Hamiltonian (general one-body h1 + up-down two-body
    Upqrs) over configurations of `nup`/`ndn` electrons in the orbital set.

    active : list of active orbital indices (default: all n). core : list of
    always-occupied orbitals (frozen, added to every determinant; their one-body
    + Hartree energy folded into the returned configs via h1/ecore externally for
    the full-space gate we pass active=all, core=[]). Returns (E0, configs) with
    configs = sorted [(coef, up_tuple, dn_tuple)] in FULL orbital indices."""
    active = list(range(n)) if active is None else list(active)
    core = [] if core is None else list(core)
    nau = nup - len(core)        # active electrons (assume core doubly occupied)
    nad = ndn - len(core)
    ups = list(combinations(active, nau))
    dns = list(combinations(active, nad))
    cfgs = [(u, d) for u in ups for d in dns]
    index = {c: m for m, c in enumerate(cfgs)}
    M = len(cfgs)
    H = np.zeros((M, M))
    core_t = tuple(sorted(core))
    # full occupations include the frozen core (always occupied, both spins)
    def full(occ):
        return tuple(sorted(core_t + occ))
    # precompute connectivity: only config pairs differing by <=1 up and <=1 dn
    for m, (u, d) in enumerate(cfgs):
        fu, fd = full(u), full(d)
        for mp, (up_, dp_) in enumerate(cfgs):
            if mp < m:
                continue
            fup, fdp = full(up_), full(dp_)
            tu = _spin_trans(fup, fu)
            td = _spin_trans(fdp, fd)
            if tu is None or td is None:
                continue
            val = 0.0
            # one-body: up changes / dn unchanged, or dn changes / up unchanged
            if td[0] == "diag":
                val += _one_body_me(tu, h1)
            if tu[0] == "diag":
                val += _one_body_me(td, h1)
            val += _two_body_me(tu, td, Upqrs)
            if m == mp:
                val += ecore
            H[mp, m] += val
            if mp != m:
                H[m, mp] += val
    w, V = np.linalg.eigh(H)
    c0 = V[:, 0]
    out = sorted(((float(c0[m]), full(cfgs[m][0]), full(cfgs[m][1])) for m in range(M)),
                 key=lambda x: -abs(x[0]))
    return float(w[0]), out


def build_multidet_trial(W, configs, ncfg, n, nup, ndn):
    """Map the top-`ncfg` CAS configs (orbital-index occupations) to site-basis
    Slater matrices for TrialWF.set_multidet: each determinant's columns are the
    occupied natural-orbital vectors W[:, occ]."""
    sub = configs[:ncfg]
    coef = np.array([c for c, _, _ in sub])
    ups = [W[:, list(u)].copy() for _, u, _ in sub]
    dns = [W[:, list(d)].copy() for _, _, d in sub]
    coef = coef / np.linalg.norm(coef)
    return ups, dns, coef


def natural_orbitals(rho):
    """Natural orbitals of a (symmetrized) 1-RDM, occupations descending."""
    occ, vecs = np.linalg.eigh(0.5 * (rho + rho.T))
    order = np.argsort(-occ)
    return occ[order], vecs[:, order]


def ensemble_rdm(cpmc):
    """Trial-agnostic charge + spin 1-RDMs from the walker ensemble: each walker
    contributes its own projector P = phi (phi^T phi)^{-1} phi^T (P_ij = <c^+_i c_j>,
    symmetric), weight-averaged. Returns (rho_up, rho_dn)."""
    wk = cpmc.walkers; n = cpmc.n
    ru = np.zeros((n, n)); rd = np.zeros((n, n)); W = 0.0
    for i in range(wk.nw):
        if wk.w[i] <= 0:
            continue
        pu = wk.phi_up[i]; pd = wk.phi_dn[i]
        try:
            gu = pu @ np.linalg.solve(pu.T @ pu, pu.T)
            gd = pd @ np.linalg.solve(pd.T @ pd, pd.T)
        except np.linalg.LinAlgError:
            continue
        w = wk.w[i]; ru += w * gu; rd += w * gd; W += w
    if W <= 0:
        return np.eye(n), np.eye(n)
    return ru / W, rd / W


def occupation_entropy(occ):
    """Sum of per-orbital entanglement entropies of the charge natural occupations
    (occ in [0,2]): S = -sum_p [f ln f + (1-f) ln(1-f)], f = occ_p/2. A scalar
    'multireference-ness': ~0 when every orbital is empty/doubly-occupied (single
    reference), large when many orbitals are fractionally occupied (the state is
    hard to converge -> a bigger active space is needed)."""
    f = np.clip(np.asarray(occ) / 2.0, 1e-12, 1 - 1e-12)
    return float(-(f * np.log(f) + (1 - f) * np.log(1 - f)).sum())


def auto_active_space(occ, nup, ndn, thr=0.05, max_active=8, min_active=2):
    """Pick (n_core, n_active) from the descending charge occupations: freeze
    orbitals with occ > 2-thr (doubly occupied), drop orbitals with occ < thr
    (empty), correlate the fractional middle. The active size GROWS automatically
    where more orbitals are fractionally occupied (strong correlation / near a
    transition / sign-problem region). Capped at max_active for cost; the cap
    being hit is itself the 'method is stretched here' signal."""
    occ = np.asarray(occ)
    n = len(occ)
    n_core = int((occ > 2 - thr).sum())
    n_virt = int((occ < thr).sum())
    n_active = n - n_core - n_virt
    n_core = min(n_core, nup, ndn)                 # need active electrons >= 0
    n_active = max(min_active, min(n_active, max_active, n - n_core))
    # keep the active electron count in range [0, 2*n_active]
    nau = nup - n_core; nad = ndn - n_core
    while (nau > n_active or nad > n_active) and n_active < n - n_core:
        n_active += 1
    return n_core, min(n_active, n - n_core)


def cas_accuracy_estimate(K, U, W_no, nup, ndn, n_core, n_active):
    """ED-free truncation-error proxy: |E_cas(n_active+1) - E_cas(n_active)|. When
    small the active space is converged (accurate); when large the wavefunction is
    multireference and the (capped) CAS is missing correlation -- flagging where the
    method is least trustworthy. Returns (E0_cas, dE_trunc, capture)."""
    n = K.shape[0]
    E0, configs = casci_solve(*transform_hubbard(K, U, W_no), n, nup, ndn,
                              active=list(range(n_core, n_core + n_active)),
                              core=list(range(n_core)))
    capture = float(sum(c * c for c, _, _ in configs))   # normalized -> ~1
    if n_core + n_active < n:
        E1, _ = casci_solve(*transform_hubbard(K, U, W_no), n, nup, ndn,
                            active=list(range(n_core, n_core + n_active + 1)),
                            core=list(range(n_core)))
        dE = abs(E1 - E0)
    else:
        dE = 0.0                                          # full active = exact
    return E0, dE, capture


def casci_trial(K, U, W_no, nup, ndn, n_core, n_active):
    """Build a CASCI multi-determinant trial in the natural-orbital basis W_no
    (columns sorted by descending occupation): freeze the top `n_core` orbitals
    (doubly occupied), correlate the next `n_active` (CAS), drop the rest. Returns
    (E0_cas, ups, dns, coef, ncfg) with site-basis Slater determinants."""
    n = K.shape[0]
    h1, Upqrs = transform_hubbard(K, U, W_no)
    core = list(range(n_core))
    active = list(range(n_core, n_core + n_active))
    E0, configs = casci_solve(h1, Upqrs, n, nup, ndn, active=active, core=core)
    ncfg = len(configs)
    ups, dns, coef = build_multidet_trial(W_no, configs, ncfg, n, nup, ndn)
    return E0, ups, dns, coef, ncfg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=2); ap.add_argument("--ly", type=int, default=2)
    ap.add_argument("--nup", type=int, default=2); ap.add_argument("--ndn", type=int, default=2)
    ap.add_argument("--U", type=float, default=4.0); ap.add_argument("--t", type=float, default=1.0)
    ap.add_argument("--basis", choices=["site", "Keig", "random"], default="Keig",
                    help="orbital basis W to test basis-invariance of CASCI E0")
    a = ap.parse_args()
    n = a.lx * a.ly
    from cpqmc import square_hopping
    K = square_hopping(a.lx, a.ly, a.t)

    if a.basis == "site":
        W = np.eye(n)
    elif a.basis == "random":
        Q, _ = np.linalg.qr(np.random.default_rng(0).standard_normal((n, n))); W = Q
    else:
        _, W = np.linalg.eigh(K)            # natural-orbital-like (eigenbasis of K)

    h1, Upqrs = transform_hubbard(K, a.U, W)
    E0, configs = casci_solve(h1, Upqrs, n, a.nup, a.ndn)
    msg = f"CASCI E0 (basis={a.basis}) = {E0:.8f}  over {len(configs)} configs"
    try:
        from hubbard_ed import build
        _, H, _, _ = build(a.lx, a.ly, a.nup, a.ndn, a.t, a.U)
        E0q = float(np.linalg.eigvalsh(H.toarray())[0])
        msg += f"   | QuSpin ED {E0q:.8f}  diff {abs(E0 - E0q):.2e}"
    except Exception as e:
        msg += f"   | (no ED cross-check: {e})"
    print(msg)


if __name__ == "__main__":
    main()
