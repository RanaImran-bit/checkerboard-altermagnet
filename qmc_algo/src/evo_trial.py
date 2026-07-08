#!/usr/bin/env python3
"""Algorithm A(a) — Evolutionary (genetic-algorithm) selected-CI trial.

Design: ../notes/00_design.md. This is the first, ED-gateable step of the
evolutionary/neural trial-state programme: GENERATE a good multideterminant
trial by EVOLUTION (selection + crossover + mutation) instead of by ED, then
hand it to the validated CPMC constraint via TrialWF.set_multidet.

Genotype:   a configuration  cfg = (up_tuple, dn_tuple)  of occupied ORBITAL
            indices in a fixed orthonormal basis W (we use the eigenbasis of K
            / natural orbitals, so the non-interacting reference is a SINGLE
            determinant and excitations are compact).
Population:  an archive (set) of configs = a selected-CI expansion. Its energy
            is the lowest eigenvalue of the selected-CI Hamiltonian over the
            archive (Slater-Condon, reusing casci.py's validated matrix-element
            helpers).
Fitness:     marginal energy lowering of a candidate when added to the archive,
            estimated by Epstein-Nesbet 2nd-order PT (the CIPSI selection
            criterion) -- cheap, NO re-diagonalization per candidate, NO ED.
Operators:   mutation = particle-hole excitation; crossover = spin-sector
            recombination + orbital swap between two parent configs.

Self-validating gates (see validate_evo_trial.py):
  (1) EXACTNESS: archive -> all configs reproduces ED (basis-invariant CI).
  (2) BIAS: feeding the evolved archive as the CPMC trial drives the
      constrained-path bias down monotonically with archive size, with NO ED
      used to build the trial -- the thesis of the design note.
"""
from __future__ import annotations
import os, sys
from itertools import combinations
import numpy as np

# reuse the validated CI matrix-element helpers + trial builder from pyqmc/casci
_HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(_HERE, "..", "..", "pyqmc"))
sys.path.insert(0, os.path.join(_HERE, "..", "..", "ed"))
from casci import (transform_hubbard, _spin_trans, _one_body_me,  # noqa: E402
                   _two_body_me, build_multidet_trial)


# --------------------------------------------------------------------------
# Selected-CI solver over an ARBITRARY set of configs (cf. casci.casci_solve,
# which enumerates a full active-space product). Same Slater-Condon convention.
# --------------------------------------------------------------------------
def ci_matrix_element(cfg_bra, cfg_ket, h1, Upqrs):
    """<bra|H|ket> for two configs (each (up_tuple, dn_tuple)) via Slater-Condon,
    general one-body h1 + up-down two-body Upqrs. Returns 0 if they differ by
    more than a single excitation in either spin."""
    ub, db = cfg_bra
    uk, dk = cfg_ket
    tu = _spin_trans(ub, uk)
    td = _spin_trans(db, dk)
    if tu is None or td is None:
        return 0.0
    val = 0.0
    if td[0] == "diag":
        val += _one_body_me(tu, h1)
    if tu[0] == "diag":
        val += _one_body_me(td, h1)
    val += _two_body_me(tu, td, Upqrs)
    return val


def selected_ci_solve(configs, h1, Upqrs):
    """Diagonalize H restricted to the given list of configs. Returns
    (E0, civec, configs) with civec the ground-state CI vector (aligned to
    `configs` order)."""
    M = len(configs)
    H = np.zeros((M, M))
    for m in range(M):
        for mp in range(m, M):
            val = ci_matrix_element(configs[mp], configs[m], h1, Upqrs)
            if val == 0.0:
                continue
            H[mp, m] += val
            if mp != m:
                H[m, mp] += val
    w, V = np.linalg.eigh(H)
    return float(w[0]), V[:, 0], configs


# --------------------------------------------------------------------------
# Fitness: Epstein-Nesbet 2nd-order energy lowering of a candidate config
# relative to the current archive ground state (the CIPSI selection criterion).
#     dE_c = |<c|H|Psi_A>|^2 / (E_A - <c|H|c>)   (<= 0, more negative = better)
# --------------------------------------------------------------------------
def en_score(cand, archive, civec, E_A, h1, Upqrs):
    """Importance of a candidate = MAGNITUDE of its Epstein-Nesbet energy
    correction, returned NEGATIVE so the most important sorts first. We take the
    magnitude (not the signed correction) deliberately: in a multireference state
    the dominant external configs can have diagonal energy BELOW the current
    archive energy (the HF reference is a poor, high-lying start), which flips the
    EN denominator sign. A sign-based filter would then reject exactly the configs
    we most need. Regularized denominator avoids a blow-up at Hcc ~ E_A."""
    coupling = 0.0
    for m, cm in enumerate(archive):
        coupling += civec[m] * ci_matrix_element(cand, cm, h1, Upqrs)
    if coupling == 0.0:
        return 0.0                                  # not (yet) connected to archive
    Hcc = ci_matrix_element(cand, cand, h1, Upqrs)
    denom = E_A - Hcc
    denom = denom if abs(denom) > 0.5 else (0.5 if denom >= 0 else -0.5)
    return -abs(coupling * coupling / denom)        # most negative = most important


# --------------------------------------------------------------------------
# Genetic operators on configs (occupation strings of orbital indices)
# --------------------------------------------------------------------------
def _excite_spin(occ, n, rng, level=1):
    """Apply `level` particle-hole excitations to one spin's occupation tuple."""
    occ = set(occ)
    for _ in range(level):
        virt = [o for o in range(n) if o not in occ]
        if not occ or not virt:
            break
        out = int(rng.choice(sorted(occ)))
        into = int(rng.choice(virt))
        occ.discard(out)
        occ.add(into)
    return tuple(sorted(occ))


def mutate(cfg, n, rng, double_prob=0.3):
    """Particle-hole excitation of a random spin sector (single, or double with
    probability double_prob)."""
    u, d = cfg
    level = 2 if rng.random() < double_prob else 1
    if rng.random() < 0.5:
        return (_excite_spin(u, n, rng, level), d)
    return (u, _excite_spin(d, n, rng, level))


def _single_excitations(cfg, n):
    """All single particle-hole excitations of a config (both spin sectors). This
    is the deterministic LOCAL-search part of the memetic algorithm: the GA's
    random crossover/mutation provides global moves, but local completeness (every
    connected single) is what makes selected-CI competitive with CIPSI/heat-bath."""
    u, d = cfg
    out = []
    for (occ, other, spin) in ((u, d, 0), (d, u, 1)):
        socc = set(occ)
        virt = [o for o in range(n) if o not in socc]
        for a in occ:
            for v in virt:
                ne = tuple(sorted(socc - {a} | {v}))
                out.append((ne, other) if spin == 0 else (other, ne))
    return out


def crossover(c1, c2, n, rng):
    """Recombine two parent configs. With prob 1/2: spin-sector recombination
    (offspring = up of one parent, dn of the other) -- always a valid config.
    Otherwise: swap one occupied orbital between the parents' up strings (kept
    valid by enforcing the particle number)."""
    u1, d1 = c1
    u2, d2 = c2
    if rng.random() < 0.5:
        return (u1, d2)
    # orbital-swap on the up sector
    s1 = set(u1)
    pool = sorted(set(u2) - s1)
    if not pool:
        return (u1, d2)
    add = int(rng.choice(pool))
    drop = int(rng.choice(sorted(s1)))
    s1.discard(drop)
    s1.add(add)
    return (tuple(sorted(s1)), d1)


# --------------------------------------------------------------------------
# The evolutionary loop
# --------------------------------------------------------------------------
def evolve_trial(h1, Upqrs, n, nup, ndn, k_target, rng,
                 ngen=40, n_candidates=64, cross_prob=0.5, seed_configs=None,
                 patience=6, local_top=6, verbose=False):
    """Evolve a selected-CI archive of up to k_target configs that minimizes the
    selected-CI energy, using ONLY the model Hamiltonian (no ED).

    Returns (E0, archive, civec): archive is a list of (up_tuple, dn_tuple),
    civec the aligned CI coefficients (descending |c| not guaranteed; caller
    sorts when building the trial)."""
    # Evolve to a WORKING size larger than k_target, then truncate to k_target at
    # the end. Growing straight to a small k_target with per-generation pruning is
    # greedy: a config that is small now but becomes important once its partners
    # arrive gets culled prematurely (observed: direct k=8 stalled at the 2-config
    # energy). The overshoot lets such configs survive; the final truncation keeps
    # the best k_target by |coefficient|.
    k_work = min(_space_size(n, nup, ndn), max(k_target + 24, 4 * k_target))
    ref = (tuple(range(nup)), tuple(range(ndn)))    # non-interacting reference
    archive = [ref]
    seen = {ref}
    if seed_configs:
        for c in seed_configs:
            c = (tuple(sorted(c[0])), tuple(sorted(c[1])))
            if c not in seen:
                archive.append(c); seen.add(c)
    E, civec, _ = selected_ci_solve(archive, h1, Upqrs)
    best_E = E
    stall = 0
    for gen in range(ngen):
        if len(archive) >= k_work:
            break
        # --- reproduce: candidate pool = LOCAL single excitations of the most
        #     important archive members (deterministic completeness) + GA offspring
        #     by crossover/mutation (global moves). ---
        cands = set()
        members = archive
        top = sorted(range(len(archive)), key=lambda i: -abs(civec[i]))[:local_top]
        for i in top:
            for c in _single_excitations(archive[i], n):
                if c not in seen:
                    cands.add(c)
        # GA offspring (bounded attempts; collisions don't advance in small spaces)
        max_tries = 30 * n_candidates
        for _try in range(max_tries):
            if len(cands) >= n_candidates + len(top) * nup * (n - nup):
                break
            if len(members) >= 2 and rng.random() < cross_prob:
                a, b = (int(i) for i in rng.choice(len(members), size=2, replace=False))
                c = crossover(members[a], members[b], n, rng)
            else:
                a = int(rng.choice(len(members)))
                c = mutate(members[a], n, rng)
            c = (tuple(sorted(c[0])), tuple(sorted(c[1])))
            if len(c[0]) == nup and len(c[1]) == ndn and c not in seen:
                cands.add(c)
        if not cands:
            break
        # --- score by |Epstein-Nesbet correction| (most negative = best) ---
        scored = sorted(((en_score(c, archive, civec, E, h1, Upqrs), c) for c in cands))
        room = k_work - len(archive)
        admitted = [c for s, c in scored[:room] if s < 0.0]
        # --- immigration: if perturbative candidates are exhausted (all
        #     uncoupled) but we still want more, inject random unseen configs so
        #     the search can escape a plateau / reach the full space (exactness). ---
        if not admitted:
            fresh = [c for _, c in scored if c not in seen][:room]
            admitted = fresh
        if not admitted:
            break
        for c in admitted:
            archive.append(c); seen.add(c)
        E, civec, _ = selected_ci_solve(archive, h1, Upqrs)
        # --- cull to the WORKING size by |CI coefficient| ---
        if len(archive) > k_work:
            order = sorted(range(len(archive)), key=lambda i: -abs(civec[i]))[:k_work]
            order = sorted(order)
            archive = [archive[i] for i in order]
            E, civec, _ = selected_ci_solve(archive, h1, Upqrs)
        if verbose:
            print(f"  gen {gen:>3}: |archive|={len(archive):>4}  E={E:.6f}")
        if E < best_E - 1e-9:
            best_E = E; stall = 0
        else:
            stall += 1
            if stall >= patience:
                break
    # --- final truncation to k_target by |CI coefficient| ---
    if len(archive) > k_target:
        order = sorted(range(len(archive)), key=lambda i: -abs(civec[i]))[:k_target]
        archive = [archive[i] for i in sorted(order)]
        E, civec, _ = selected_ci_solve(archive, h1, Upqrs)
    return E, archive, civec


def _space_size(n, nup, ndn):
    from math import comb
    return comb(n, nup) * comb(n, ndn)


# --------------------------------------------------------------------------
# Build a CPMC-ready multideterminant trial from the evolved archive
# --------------------------------------------------------------------------
def evolved_multidet(W, archive, civec):
    """Map the evolved archive (configs in orbital basis W) + CI coefficients to
    site-basis Slater determinants for TrialWF.set_multidet, sorted by |coef|."""
    sub = sorted(((float(civec[m]), archive[m][0], archive[m][1])
                  for m in range(len(archive))), key=lambda x: -abs(x[0]))
    return build_multidet_trial(W, sub, len(sub), W.shape[0], None, None)
