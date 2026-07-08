# Equation reference — every formula used by the platform's engines & estimators

One place for the math implemented across ED / DQMC / CP-DQMC / CPQMC and the
chi_zz(q) + AHE work. Each block cites the implementing code. Conventions are the
CODE's conventions (they differ between engines — noted explicitly).

---

## 1. Model (`pyqmc/cpqmc.py::am_hopping`, `square_hopping`; K = −hopping)

    H = Σ_σ Σ_ij K_σ[i,j] c†_iσ c_jσ + U Σ_i n_i↑ n_i↓  (− μ Σ_i n_i at finite T)

Altermagnet dispersion (up spin; down = 90°-rotated, x↔y / diagonals swapped):

    ε_↑(k) = −2(t0−tam)cos kx − 2(t0+tam)cos ky  − 2 t1[cos(kx+ky) − cos(kx−ky)] − 2 tp[cos(kx+ky)+cos(kx−ky)]
    ε_↑(k) − ε_↓(k) = 4 tam (cos kx − cos ky) − 4 t1 [cos(kx+ky) − cos(kx−ky)]      (d_{x²−y²} + d_{xy} spin splitting)

Site index everywhere: `i = x*ly + y` (PBC). Two-orbital variant: `i = orb*lxy + x*ly + y`.

## 2. Trotter + Hubbard-Stratonovich

Finite T (`code/dqmc_py/dqmc.py`, `code/ftcpqmc_py/ftcpmc.py`): β = L·Δτ,

    e^{−βH} = Π_{l=1..L} e^{−Δτ K'} e^{−Δτ V_l} + O(Δτ²),   μ_eff = μ − U/2 folded into K'

Discrete Hirsch HS (particle-hole symmetric, x_{i,l} = ±1):

    e^{−Δτ U (n↑−½)(n↓−½)} = ½ Σ_{x=±1} e^{λ x (n↑ − n↓)} · e^{−Δτ U/4},   cosh λ = e^{Δτ U / 2}

Per-slice propagator:

    B_l^σ = e^{−Δτ (K_σ − μ_eff)} · diag(e^{s_σ λ x_{i,l}}),   s_↑ = +1, s_↓ = −1

T=0 CPQMC (`pyqmc/cpqmc.py::LatticeModel/Propagator`): symmetric split
b^σ = e^{−ΔτK/2} D^σ(x) e^{−ΔτK/2}; general density-density terms (uxy, v) via the
same discrete HS in difference/sum channels (`LatticeModel._hs`, `dqmc._dd_terms`):

    V>0: e^{−ΔτV n_a n_b} = e^{−ΔτV(n_a+n_b)/2} ½Σ_s e^{λ s (n_a − n_b)},  cosh λ = e^{ΔτV/2}
    V<0: same with (n_a + n_b) and cosh λ = e^{−ΔτV/2}

## 3. Finite-T determinant QMC (BSS) (`code/dqmc_py/dqmc.py`)

    Z = Σ_{x} Π_σ det( I + B_L^σ B_{L−1}^σ … B_1^σ )
    G^σ(l,l) = [ I + B_l…B_1 B_L…B_{l+1} ]^{−1}       (equal-time at slice l; G_ij = <c_i c†_j>)

Metropolis single-site flip ratio + Sherman–Morrison update (Δv = e^{−2 s_σ λ x} − 1):

    R_σ = 1 + (1 − G^σ_ii) Δv,   R = R_↑ R_↓ ;   G' = G − (e_i − G_{:,i}) G_{i,:} · Δv / R_σ

Time-displaced Green's functions (`_green_tau`, `spin_susc._green_0tau`):

    G(τ_l, 0) =  [ B(l,0)^{−1} + B(L,l) ]^{−1}        = <c(τ_l) c†(0)>
    G(0, τ_l) = −[ B(L,l)^{−1} + B(l,0) ]^{−1}        = −<c†(τ_l) c(0)>ᵀ
    slice propagation: G(l+1,0) = B_l G(l,0);  G(0,l+1) = G(0,l) B_l^{−1};  G(l+1,l+1) = B_l G(l,l) B_l^{−1}

UDV / ASvQRD stabilization (`_udv_chain`; docs/STABILIZATION.md), big/small split
D = D_b·D_s, |D_b| ≥ 1 ≥ |D_s|:

    G = V2^{−1} [ U1ᵀ V2^{−1} / D1_b D2_b + D1_s (V1 U2) D2_s ]^{−1} D1_s V1 / D2_b   (schematic, code exact)

Average sign: <s> = Σ sgn(det·det) / Σ 1.

## 4. T=0 constrained-path AFQMC (CPQMC) (`pyqmc/cpqmc.py`)

    |ψ0⟩ ∝ lim_{n→∞} (e^{−ΔτH})ⁿ |ψT⟩,   walkers: |Φ⟩ = Slater det., weight w

Overlap & single-particle Green (per spin; code returns g = <c†_i c_j>):

    O = det(ΨTᵀ Φ),   g = [ Φ (ΨTᵀ Φ)^{−1} ΨTᵀ ]ᵀ

Importance sampling + constrained path: field chosen from heat-bath of overlap
ratios; walker killed when ⟨ψT|Φ⟩ ≤ 0. Mixed and back-propagated estimators:

    <O>_mix = ⟨ψT| O |Φ⟩ / ⟨ψT|Φ⟩ ;    <O>_BP = ⟨ψT B_bp…B_1 | O | Φ(0)⟩ / ⟨·|·⟩

(BP bra from `bp_bra` over the recorded fields; kets frozen at the block start.)

Time-displaced GFs inside the BP window (`chi_block`/`chi_spin_block`/`chid_block`),
with B_(l) = b_{l}…b_1 the recorded one-body propagators and B^{−1} accumulated:

    particle  P^σ(τ_l) = B_(l) (I − g^{σT})   = <c_i(τ_l) c†_j(0)>
    hole      H^σ(τ_l) = B_(l)^{−T} g^σ       = <c†_i(τ_l) c_j(0)>
    equal-time g^σ(τ_l) = B^{−T} g^σ B^{T}    (diagonal → n_i(τ_l))

## 5. Finite-T constrained-path AFQMC (CP-DQMC) (`code/ftcpqmc_py/ftcpmc.py`)

Trial-propagator importance function (b0 = U=0 propagator, T_l = b0^{L−l}):

    I_l = Π_σ det( I + T_l^σ M_l^σ ),   M_l^σ = B_l^σ…B_1^σ
    I_0 = Z_T,  I_L = W(x)  ⇒  Π_l r_l = W(x)/Z_T,  r_l = I_l / I_{l−1}

Constraint (Zhang finite-T CP): kill walker when a branch ratio ≤ 0; released = exact.
Force-biased per-site heat-bath with the local ratio (same R as DQMC), weight
log w += log[(|f₊|+|f₋|)/2]; stabilized Green's from the DQMC ASvQRD chain.

## 6. Wick contractions / equal-time estimators

CPQMC convention g_ij = <c†_i c_j> (n_i = g_ii); DQMC convention G_ij = <c_i c†_j>
(n_i = 1 − G_ii, <c†_i c_j> = δ_ij − G_ji). Same-spin exchange from one HS config:

    <n_i n_j>   = (n↑+n↓)_i (n↑+n↓)_j + Σ_σ [δ_ij n^σ_i − G^σ_ij G^σ_ji]
    <S^z_i S^z_j> = ¼ [ (n↑−n↓)_i (n↑−n↓)_j + Σ_σ (δ_ij n^σ_i − G^σ_ij G^σ_ji) ]
    energy: E = Σ_σ Σ_ij K_σ[i,j] <c†_i c_j>_σ + U Σ_i n↑_i n↓_i (+ uxy/v Wick terms)

Singlet pairing (Δ†_α = Σ_m Σ_δ f_α(δ) c†_{m↑} c†_{m+δ↓}; F_α bond form-factor matrix,
f_s = (+1,+1,+1,+1), f_d = (+1,+1,−1,−1) on (±x, ±y)):

    P_α(m,n) = <Δ_α(m) Δ†_α(n)> = G↑[m,n] · (F_α G↓ F_αᵀ)[m,n]     (FULL)
    P_α^vertex = P_α − Ḡ↑ ∘ (F_α Ḡ↓ F_αᵀ)                          (bubble from ensemble-avg G)

## 7. Dynamic (imaginary-time) susceptibilities

Pairing (`chid_block`, `dqmc.chi_pair`, `ftcpmc chi=True`):

    C_α(τ) = Σ_{m,n} P^↑(τ)_{mn} (F_α P^↓(τ) F_αᵀ)_{mn},   χ_α = ∫_0^{β or τmax} dτ C_α(τ)
    χ_α^vertex = χ_α − ∫dτ Σ P̄^↑(τ)∘(F_α P̄^↓(τ)F_αᵀ)      (τ-resolved ensemble-avg bubble)

**Spin (the new chi_zz(q); `spin_susc.py`, `ftcpmc spin=True`, `cpqmc.chi_spin_block`):**
per-configuration Wick (cross-spin factorizes; disconnected term SURVIVES config
averaging and carries the q-structure):

    <S^z_i(τ) S^z_j(0)>_cfg = ¼ [ m_i(τ) m_j(0) − Σ_σ G_σ(0,τ)_{ji} G_σ(τ,0)_{ij} ],  m_i = n_{i↑} − n_{i↓}

CPQMC (T=0) matrix form (identical content in the <c†c> convention):

    M_l[i,j] = m̃_i(τ_l) m̃_j(0) + ¼ Σ_σ (H^σ(τ_l) ∘ P^σ(τ_l))[i,j],   m̃ = ½(n↑ − n↓)

Momentum reduction (`reduce_mat` / `_shift_index`; identical across engines):

    S(R) = Σ_m M[m, m+R],   P(q) = Re FFT2[S(R)],   χ_zz(q) = P(q)/N,  q = 2π(kx/lx, ky/ly)

Integration: finite-T rectangle rule Δτ Σ_l (matches chi_pair); T=0 windowed trapezoid

    χ^win_q = Δτ [ ½C_q(0) + Σ_{0<l<bp} C_q(τ_l) + ½C_q(τ_bp) ]

**Convention bridge (ED-verified):**  χ^{finite-T}(β→∞) = 2 · χ^{T=0, one-sided}
(∫_0^β picks up both time orderings);  equal-time S^z(q) identical in both.

Pairing eigenvalue (`pair_eig.py`, `ftcpmc paireig_tau`): connected k-space pair matrix

    P_c(k,k′) = ∫dτ [ <G↑(k,k′;τ) G↓(−k,−k′;τ)> − <G↑><G↓> ],  leading eig λ, eigvec φ(k)
    d/s character: |⟨f_{d,s}|φ⟩|,  f_d = cos kx − cos ky, f_s = cos kx + cos ky

## 8. Exact diagonalization references

Sector ED (QuSpin `hubbard_ed.build` / numpy `validate_chi_spin._sector_hop`):
fixed (N↑, N↓) bitmask basis, Jordan-Wigner sign for c†_i c_j = (−1)^{#occupied strictly between i,j};
H = H↑ ⊗ I + I ⊗ H↓ + U Σ_i n↑_i n↓_i. Full-Fock JW (`spin_susc._fock_ops`): 2^{2N} space,
c_a = parity(occ_{<a}) on bit a.

Finite-T Lehmann/Kubo (any operator O_q; `spin_susc.ed_chi_spin`, `dqmc.ed_finite_T`):

    χ_O(q) = (1/ZN) Σ_{a,b} |⟨b|O_q|a⟩|² K(E_a,E_b),
    K(E_a,E_b) = (e^{−βE_b} − e^{−βE_a}) / (E_a − E_b),   K → β e^{−βE_a} as E_b → E_a
    S_O(q)  = (1/ZN) Σ_a e^{−βE_a} ⟨a|O_q O_q†|a⟩,       O_q = Σ_i e^{−iq·r_i} S^z_i

T=0 Lehmann (`validate_chi.ed_chi`, `validate_chi_spin.ed_chi_spin_T0`):

    C_q(τ) = (1/N) Σ_n |⟨n|O_q|0⟩|² e^{−(E_n−E_0)τ},   χ_q = (1/N) Σ_{n>0} |⟨n|O_q|0⟩|²/(E_n−E_0)

## 9. Statistics

Block averages: value = mean over blocks/paths; error = σ_blocks/√N_blocks.
Sign-weighted averages: <O> = Σ w s O / Σ w s, <sign> = Σ s / Σ |1|.
Regression gate: fixed seed ⇒ identical RNG stream ⇒ bit-for-bit equality required.

## 10. AHE formulas (strategy stage; docs/AHE_MEASUREMENT.md)

    Currents:      j_x = i Σ_{⟨ij⟩_x,σ} K_σ[i,j] (c†_i c_j − c†_j c_i);   j_x^z = same weighted by s_σ
    Kubo:          σ_xy = lim_{ω→0} (1/ω) Im Λ^A_{xy}(ω),  Λ_{xy}(iω_n) = ∫_0^β dτ e^{iω_n τ} ⟨j_x(τ) j_y(0)⟩,
                   Λ^A = ½(Λ_xy − Λ_yx);   first-Matsubara proxy σ_xy ≈ Λ^A(iω_1)/ω_1
    Spin-splitter: σ^z_xy from ⟨j_x^z(τ) j_y(0)⟩ (allowed WITHOUT SOC; charge σ_xy ≡ 0 for real K_σ)
    Streda:        σ_xy = e² ∂n/∂B |_μ  (Peierls: K_ij → K_ij e^{i(e/ħ)∫A·dl}, Φ = 2π/N_φ)
    Top. Ham.:     h_top(k) = −G(k, iω→0)^{−1};  Chern # via lattice field strength (FHS)
    Chern marker:  C(r) = −2π Im Tr { P [x̂,P][ŷ,P] },  P = one-body density matrix
