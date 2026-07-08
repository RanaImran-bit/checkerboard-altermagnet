# PLAN — Momentum-space magnetic susceptibility χ_zz(q) in all engines + AHE strategy

**Date:** 2026-07-08 · **Scope:** magnetic measurements of the altermagnet Hubbard model
(single-band, `am_hopping(lx,ly,t0,tam,t1,tp)`).

Two asks:
1. **χ_zz(q)** — the momentum-resolved, τ-integrated (Kubo) magnetic susceptibility — in
   **all** production engines (DQMC, CP-DQMC, CPQMC), each **gated against ED** in the
   platform's usual style (exact U=0 gates + interacting small-cluster gates).
2. A method strategy for measuring the **anomalous Hall effect** with these engines
   (→ `docs/AHE_MEASUREMENT.md`).

---

## 1. Definitions & conventions (fixed up front, all engines identical)

Finite-T (DQMC, CP-DQMC), per site:

    χ_zz(q) = (1/N) Σ_ij e^{-iq·(r_i - r_j)} ∫_0^β dτ <S^z_i(τ) S^z_j(0)>,   S^z_i = (n_iu - n_id)/2
    S^z(q)  = (1/N) Σ_ij e^{-iq·(r_i - r_j)} <S^z_i S^z_j>                    (equal-time SF)

T=0 (CPQMC), windowed one-sided Kubo integral over the BP window (same convention as the
existing `run_bp_chi` / `validate_chi` staggered gate):

    C_q(τ)  = (1/N) Σ_n |<n|S^z_q|0>|^2 e^{-(E_n-E_0)τ},   S^z_q = Σ_i e^{-iq·r_i} S^z_i
    χ_q^win = trapezoid of C_q(τ) over τ ∈ [0, bp·dt]

**Cross-convention factor (verified vs ED):** χ_finite-T(β→∞) = **2 ×** χ_T=0^one-sided
(the finite-T ∫_0^β picks up both time orderings). Equal-time S^z(q) grids agree exactly.

Numerics: site index `i = x*ly + y`; q-grids produced by the shared `reduce_mat`/FFT2
reduction (`code/dqmc_py/dqmc.py`), i.e. **identical q-grid convention across engines** —
this is what already guarantees DQMC-vs-CPQMC comparability for the pairing channel.
Per-configuration Wick (same-spin exchange only; cross-spin factorizes):

    <Sz_i(τ)Sz_j(0)>_cfg = ¼ [ m_i(τ) m_j(0) − Σ_σ G_σ(0,τ)_{ji} G_σ(τ,0)_{ij} ]

The DISCONNECTED m·m term survives config averaging and carries the q-structure — never drop it.

## 2. Status quo (what already exists — do not rebuild)

| Engine | χ_zz(q) status | Reuse |
|---|---|---|
| **DQMC** | **DONE** — `code/dqmc_py/spin_susc.py::SpinDQMC` (full q-grid, stable G(τ,0)/G(0,τ)/G(l,l) with UDV restabilization) + its own full-Fock JW ED oracle `ed_chi_spin` (finite-T Lehmann/Kubo) | the Wick formula + ED oracle are the templates for the other two engines |
| **CPQMC (T=0)** | partial — `Estimators.chi_block` measures C(τ) only for a **hard-coded (π,π) staggered phase** (scalar) | its particle/hole propagators P^σ(τ)=B(I−g^T), H^σ(τ)=B^{-T}g and the `_step_bmats` machinery generalize verbatim to the full matrix |
| **CP-DQMC** | none (only equal-time S_AM) | the `chi=True` branch of `_accum` already walks the τ-slices with stable `_green_tau` restabilization — add G(0,τ) and the spin Wick term |

## 3. Work items

**W0 — git baseline.** The tree lost its `.git` (unzipped copy). `git init`, ignore the
zip, commit the existing tree as baseline so every change below is a reviewable diff.

**W1 — CPQMC (T=0) general-q spin susceptibility.** `pyqmc/cpqmc.py`:
- `Estimators.chi_spin_block(walkers, ket_up, ket_dn, rec, bp)` — exact matrix
  generalization of `chi_block`: per walker, per τ-slice l, accumulate
  `M_l[i,j] = m_i(τ_l) m_j(0) + ¼ Σ_σ (H^σ ∘ P^σ)[i,j]` (with m = ½(n_u−n_d); the phase
  contraction `phase·M·phase` of the old scalar code is deferred to the FFT). Returns
  `(Msum[L,n,n], W)`.
- `CPMC.run_bp_chi_spin(...)` — block driver mirroring `run_bp_chi`; per block reduce each
  `M_l` with `reduce_mat` → `C_q(τ_l)` grids; outputs `taus, Ctau_q (L,lx,ly) ± err`,
  windowed `chi_q (lx,ly) ± err`, plus the (π,π) trace for continuity with `run_bp_chi`.
- **Gate** `pyqmc/validate_chi_spin.py`: T=0 Lehmann ED oracle at general q (QuSpin sector
  ED, extends `validate_chi.ed_chi` to complex phases), compared q-by-q and τ-by-τ.
  Points: 2×2 half-filled U=4 (tam=0, 0.3); 4×2 U=0 exact gate; 4×2 half-filled U=4.
- **Regression:** additive methods only — `pyqmc/regression.py` must stay bit-for-bit.

**W2 — CP-DQMC spin susceptibility.** `code/ftcpqmc_py/ftcpmc.py`:
- Make `self._dq` a `SpinDQMC` (drop-in subclass; adds `_green_0tau`).
- `spin=True` path in `_accum`: port `SpinDQMC.chi_spin` verbatim — propagate
  Gl0/G0l/Gll per slice (incl. dd-field Bv factors), restabilize every `nstab` via
  `_green_tau/_green_0tau/green`, accumulate `Mzz` (τ-integrated) and `Szz` (equal-time)
  sign-weighted.
- `_finalize`: `chi_spin_q = reduce_mat(Mzz/W)["Pq"]/n`, same for `S_spin_q`; scalars
  `chi_q0`, `chi_max`. CLI `--spin` (+ `--ed` prints the `ed_chi_spin` grids).
- **Gate:** free-projection (`--free`) must reproduce ED exactly (within stats); the
  constrained run shows the CP bias and <sign>.

**W3 — three-way finite-T benchmark.** `code/dqmc_py/validate_chi_spin_ft.py`: one driver
that runs ED / SpinDQMC / FTCPMC(free) / FTCPMC(CP) on the same points and prints the
grids + max |rel. dev| + PASS/FAIL. Points: 2×2 U=4 μ=2 β=2, tam ∈ {0, 0.3}, t1 ∈ {0, 0.2};
3×2 U=4 β=2 tam=0.2 (JW full Fock 4^6 = 4096 — still exact). Tolerance: 3σ or 5% (whichever
larger) for DQMC/free; CP bias reported, not gated.

**W4 — AHE strategy** → `docs/AHE_MEASUREMENT.md` (suggestions + physics constraints; no
heavy implementation this round). Headline points: the current real-hopping altermagnet
has σ_xy ≡ 0 (needs SOC/complex hopping); the symmetry-allowed transverse response
*without* SOC is the **spin-splitter** channel; ranked QMC routes = topological-Hamiltonian
/ Chern marker from G, Streda dn/dB with Peierls flux (phaseless CP), antisymmetric
Matsubara current-current correlator (DQMC, sign-free at half-filling), spin-current
correlators (direct extension of this χ_zz machinery to bond-current vertices).

**W5 — simulation tests on the 48-thread server** (not the local Mac): deploy, build env
(numpy/scipy/quspin), run `pyqmc/regression.py` + all gates in parallel
(`OMP_NUM_THREADS=1` per process, many processes), pull outputs into
`results/chi_spin_bench/` and record the numbers in `docs/HANDOFF.md`.

## 4. Deliverables / gates summary

- Code: W1 + W2 (additive; regression-clean).
- Gates: `validate_chi_spin.py` (T=0), `validate_chi_spin_ft.py` (finite-T, three-way),
  existing `spin_susc.py --ed` (DQMC) — all run on the server, outputs committed.
- Docs: this plan, `AHE_MEASUREMENT.md`, HANDOFF update.
- Git: baseline commit + one commit per work item.

## 5. Physics we expect (sanity anchors)

- Half-filling, tam=0: χ_zz(q) peaks at q=(π,π) (AFM), grows with β.
- tam>0 preserves (π,π) nesting (NN d-wave splitting) → (π,π) peak survives; t1>0 breaks
  nesting → peak suppressed/shifted. The equal-time S^z(q) analogue of this is the
  established altermagnet story; χ_zz(q) adds the coherence (Kubo) weight.
- CP-DQMC at half-filling is sign=1 for all tam → χ_zz(q) in the interesting regime with
  no sign artifact, same argument as for the pairing channel.
