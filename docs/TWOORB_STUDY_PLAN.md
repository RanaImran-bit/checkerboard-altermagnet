# Two-orbital altermagnet: 4-method study plan (pairing + altermagnetism)

Model: two-orbital (d_xz,d_yz) Hubbard/altermagnet, `ed/altermagnet_ed.build_hopping`.
FIXED anisotropic bands (leeb2024spontaneous / manuscript Fig.S9): t1=1, t2=1.75, t3=0.85,
t4=0.65. Interaction is the scan axis. All four methods: ED, CPQMC (T=0), CP-DQMC (finite-T),
DQMC (finite-T).

## 0. Central questions
- PAIRING: which channel dominates (d_{x2-y2} vs ext-s); is it enhanced by interaction as an
  equal-time VERTEX correlation AND a tau-integrated SUSCEPTIBILITY; short-range or long-range?
- ALTERMAGNETISM: does AM order form SPONTANEOUSLY as interaction grows (threshold, per leeb)?
  does d-wave turn on WITH it? what competes (Neel, CDW)?
- Decisive lens for both: FINITE-SIZE SCALING (small clusters already gave artifacts).

## !! Key refinement (why S_AM looked flat in the first U-scan) !!
The first U-scan used ON-SITE U ONLY (uxy=v=0). leeb/the manuscript drive spontaneous AM with
an EXCHANGE interaction J (inter-orbital / bond S.S). Plain on-site U favors local moments but
NOT the orbital-differential magnetization that AM is -> flat S_AM is EXPECTED for the U-only
reduction. To test spontaneous AM properly the study MUST:
  (i)  add the inter-orbital interaction uxy (CPQMC & ED support it; DQMC/CP-DQMC do NOT),
  (ii) measure the AM SUSCEPTIBILITY chi_AM (tau-integrated) -- sensitive to incipient order,
  (iii) reach stronger coupling (J~4t^2/U is small at intermediate U),
  (iv) finite-size scale S_AM/N and chi_AM.
Method limitation: the inter-orbital-interaction AM study is CPQMC(+ED) only; the finite-T
DQMC/CP-DQMC codes implement on-site U only (adding uxy/v HS channels is a separate task).

## 1. Observables (consistent across methods)
PAIRING
- equal-time S_d, S_s : FULL and connected VERTEX pair structure factors
- pair correlation vs distance R (ODLRO / long-range test); size scaling S_d/N
- tau-integrated chi_d, chi_s (vertex) [finite-T; CPQMC via BP tau-window]
- leading pairing eigenvalue (form-factor-free) + d / ext-s overlap
ALTERMAGNETIC
- S_AM : orbital-staggered S^z structure factor  < (sum_i g_i S^z_i)^2 >, g=+1/-1 by orbital
  (per-site normalization in QMC = S_AM/N; ED reports extensive <O_AM^2> -> differ by N)
- chi_AM : tau-integrated AM susceptibility (the incipient-order probe)
- d-wave-projected spin splitting Delta n(k) = <n_{k up} - n_{k dn}> . (cos kx - cos ky)
- competing orders: Neel S(pi,pi), CDW structure factor, uniform moment^2
SCALING: S_d/N, S_AM/N, chi_d, chi_AM vs L (6,8,10,12) -> thermodynamic-limit behaviour.

## 2. Methods & roles (+ limits)
- ED       : exact ground truth; defines & validates every observable. Limit: <= 2x2 two-orbital
             (or small fillings).
- CPQMC    : PRIMARY, genuine two-orbital ground state at L=6..12; supports uxx/uxy/v.
             Limit: constraint bias on pairing (use trends; ED-anchor magnitudes).
- CP-DQMC  : finite-T, sign-controlled; tau-integrated chi. Limit: on-site U only; slow (python
             loops) -> needs cheap settings. NaN blow-up possible at large N,beta.
- DQMC     : finite-T, exact benchmark + chi. Limit: on-site U only; sign problem for t3,t4!=0
             (non-bipartite) -> sign-limited; monitor <sign>.

## 3. Scan axes (fixed bands)
1. Interaction: U (0..~12) on-site; AND for CPQMC/ED, uxy (0..~few) inter-orbital -- the AM driver.
2. System size L = 6,8,10,12  (finite-size scaling -- THE crux).
3. Temperature beta (finite-T methods) -- onset / T->0.
4. (optional) Filling -- the manuscript's d-wave ENHANCEMENT dome is at DOPING, not half-filling.

## 4. Phased execution
Phase 0  Infrastructure & validation
  - finalize AM observables: S_AM (done), add chi_AM (tau-integrated) + d-wave-projected Dn(k);
    fix S_AM normalization convention (per-site vs extensive) & document.
  - fix CP-DQMC speed (cheap settings: small beta, few walkers/measurements, large nstab).
  - ED-gate EVERY observable at small clusters (energy + pairing + AM). [half done]
Phase 1  U-scan at fixed size (6x6, 8x8)
  - on-site U scan: pairing vertex + chi + S_AM vs U (all 4 methods). [CPQMC/DQMC done;
    CP-DQMC pending cheap relaunch]
  - CPQMC/ED: ALSO scan uxy (inter-orbital) at fixed U -> test the AM driver directly.
Phase 2  Finite-size scaling (L=6,8,10,12) at key (U, uxy)
  - scale S_d/N, S_AM/N, chi -> thermodynamic-limit verdict on d-wave dominance & spontaneous AM.
    [10x10, 12x12 launched]
Phase 3  Temperature series (CP-DQMC/DQMC), fixed U,L
  - chi_d(T), chi_AM(T), onset temperatures.
Phase 4  (optional) Doping axis at fixed U
  - test the doped d-wave dome (the manuscript's actual enhancement claim).
Phase 5  Synthesis
  - figures: pairing & AM vs U (and vs uxy); vs 1/L; phase diagram. Findings doc + appeal note.

## 5. Validation strategy
ED-gate every observable at small size; require CPQMC(T=0) <-> CP-DQMC(large beta) <-> DQMC(sign~1)
agreement where regimes overlap; monitor DQMC <sign> and flag sign-limited points; report
per-block + seed error bars. Spontaneous SSB is not literal in a finite system -> detect via the
structure factor / susceptibility SCALING with L, not a single-size value.

## 6. Resources / logistics
Cluster boxes 251-255 (5 x 64 cores, 503 GB RAM; 254/255 had miniconda+numpy installed this
session). SSH intermittent -> jobs launched detached (setsid nohup), sharded, collected by polling.
ED runs locally. Code on each box under ~/qmc_twoorb; drivers pyqmc/{ed_uscan,two_orb_uscan}.py,
code/dqmc_py/two_orb_ft_scan.py; launchers scratch_launch_*.sh.

## 7. Deliverables
Per-phase CSVs in results/dqmc_scan/; finite-size figure set; docs/TWOORB_UAM_FINDINGS.md; the
honest verdict -- does interaction spontaneously drive AM + d-wave in the thermodynamic limit, and
is d-wave the leading channel?

## 8. Open decisions
1. Scope: core (Phases 0-2) vs full (add Phase 3 temperature / Phase 4 doping)?
2. Filling: half-filling only, or include the doped dome (Phase 4, where the manuscript's
   enhancement actually lives)?
3. AM driver: add uxy to the CPQMC/ED scans (recommended -- without it AM likely won't form)?
4. CP-DQMC: full finite-T cross-check, or spot-check a few (U,L) given its cost?

## Interaction support across methods (updated)
- ED     : full uxx/uxy/v (exact). OK.
- CPQMC  : full uxx/uxy/v (LatticeModel), ED-validated. OK.
- DQMC   : full uxx/uxy/v ADDED + ED-gated (finite-T grand-canonical sector sum: <H>,<n> match at
           sign~0.95). NB sign-limited for t3,t4!=0. OK (commit 6ff74bd).
- CP-DQMC: on-site U only. Adding uxy/v is NOT a simple port -- its trial-based stabilizer
           (future = non-interacting b0) needs a field-OFF reference (factor 1), which the ported
           _hs decoupling lacks. Correct fix = DIFFERENCE-channel decoupling
           e^{-dtV n_a n_b} = e^{-dtV(n_a+n_b)/2} (1/2) sum_s e^{lam s (n_a - n_b)} (s=0 -> factor 1),
           with the constant e^{-dtV(n_a+n_b)/2} folded into K (shifts filling). Real work + its own
           ED gate. DECISION PENDING: do the rework, or standardize on ED+CPQMC+DQMC (full model) and
           keep CP-DQMC as a U-only spot-check.

## Status snapshot (as of writing)
- U-scan (on-site U, fixed bands) at 6x6, 8x8: CPQMC + DQMC DONE; CP-DQMC pending cheap relaunch.
  10x10, 12x12: launched (CPQMC + DQMC).
- Early result: d-wave vertex GROWS with U (good) but ext-s vertex DOMINATES; S_AM flat (expected
  for U-only -- see key refinement above).
- Fixes queued: cheap CP-DQMC relaunch (scratch_launch_cpdqmc.sh); ED-vs-CPQMC S_AM cross-check
  at 4x4/nup=2; add uxy scan + chi_AM.
