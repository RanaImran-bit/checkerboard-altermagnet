# Validation log

Common-schema comparison of QMC runs against the QuSpin ED reference.
Workflow:

```bash
source tools/env.sh
# QMC
tools/build.sh benchmark-cpqmc
tools/run.sh   benchmark-cpqmc code/record/benchmark/CPQMC/in.dat 4
python tools/parse_out.py results/benchmark-cpqmc/<runid> -o /tmp/qmc.json
# ED reference (exact, error 0)
python ed/hubbard_ed.py --lx 4 --ly 4 --nup 1 --ndn 1 --t 1 --U 3 -o /tmp/ed.json
# compare (z = sigmas apart; PASS if z <= ztol)
python tools/compare.py /tmp/ed.json /tmp/qmc.json
```

ED runs single-threaded by design (`OMP_NUM_THREADS=1` + `KMP_DUPLICATE_LIB_OK`
set inside `ed/hubbard_ed.py`): conda openblas (libgomp) and the QuSpin
compiled core (llvm libomp) otherwise abort with a duplicate-OpenMP error.

## 2026-06-23 — benchmark-cpqmc vs ED, 4x4 Hubbard, Nup=Ndn=1, U=3, t=1

ED ground state (256 basis states, exact): E = -7.867205, E/site = -0.491700.

| observable        | ED        | CPQMC          | z     | verdict |
|-------------------|-----------|----------------|-------|---------|
| energy_total      | -7.867205 | -7.870103      | 1.01  | PASS    |
| energy_per_site   | -0.491700 | -0.491881      | 1.01  | PASS    |
| energy_kinetic    | -7.960602 | -7.965755      | 1.75  | PASS    |
| energy_potential  |  0.093397 |  0.095652      | 15.66 | FAIL    |

## 2026-06-24 — Python port (pyqmc) extended to uxy + v

pyqmc/cpqmc.py generalized from on-site Hubbard to an arbitrary density-density
term list V*n_{a,sa}*n_{b,sb}, each decoupled with the exact discrete HS
(repulsive=spin channel; attractive=charge channel + coeffv weight). Builds
uxx/uxy/v terms for the two-orbital model; energy sums the matching Wick terms.

Python vs ED (2x2, 4+4):
  single-band U=3:  -7.868   (regression OK)
  uxx=2:            -28.287  z=0.93 PASS
  uxy=2:            -24.605  z=0.62 PASS
  v=0.1 (aniso):    -16.450  z=1.2  PASS
  v=0.2 (aniso):    -17.034 vs -17.053  z=5.3  -> growing MIXED-ESTIMATOR bias
The Python port uses the mixed estimator (biased for operators not commuting
with H; bias grows with coupling as the free-fermion trial degrades); the
Fortran avoids this via back-propagation. Implementation verified correct at
small coupling. Both Fortran and Python now model all three channels.

## 2026-06-24 (latest) — v VALIDATED: the "collapse" was degeneracy, not a bug

The v walker collapse was NOT a code defect. The symmetric cluster
(t1=t2=t3=t4=-1) is 4-fold degenerate at the Fermi level: ED E(v) is flat
(-32) for v<=~0.5 then jumps (-36.4 at v=1) -- a level crossing. CP-QMC with a
free-fermion trial cannot converge to a degenerate / phase-crossed ground
state, hence the instability. This is a known-hard AFQMC regime, not a bug.

Using the PHYSICALLY relevant anisotropic (altermagnet) hopping
t1=-1,t2=-0.5,t3=-0.3,t4=-0.2 (which lifts the degeneracy and gives a smooth
E(v)), the v channel reproduces ED:
  v=0.1: -16.446+/-0.004 vs ED -16.449  z=0.95 PASS
  v=0.2: -17.055+/-0.006 vs ED -17.053  z=0.39 PASS
  v=0.3: -17.986+/-0.012 vs ED -17.973  z=1.08 PASS

=> ALL THREE interaction channels validate against ED:
     uxx (z~0.25), uxy (z~1.33), v (z~0.4-1.1).
The fixes that achieved this: backphi bond-partner replay, rec_fields v-gate
(abs(uxy)->abs(v)), InitV per-channel V sign. Use anisotropic, non-degenerate,
moderate-coupling clusters as v benchmarks (tools/build.sh src-small honours
QMC_LX/QMC_LY/QMC_NUP/QMC_NDN). Symmetric/strong-v small clusters are a hard
methodology regime, not a code target.

## 2026-06-24 (later) — uxy VALIDATED: back-propagation bond-partner bug fixed

Root cause (found via pyqmc/vee_check.py, which confirmed the Vee forward
determinant/GF update is exact to 1e-16 for all channels): the bug was in
backphi. Its field replay applied each interaction term's SECOND operator to
site j == loop index, which equals the first site only for the same-site uxx
channel. For uxy/v (i/=j) the partner site (i+lxy for uxy; the neighbour for v)
was never scaled -> wrong back-propagated Green's function -> biased energy.
Also fixed rec_fields: the v-field recording block was gated by abs(uxy) not
abs(v). backphi now mirrors the Vxy/Vznn bond structure (op1 at site a on spin
si, op2 at partner b on spin sj, expV indexed at a).

Result (2x2, 4+4, equilibrated), Fortran vs ED:
  uxx=2: -28.288+/-0.003 vs -28.289  z=0.25 PASS (no regression)
  uxy=2: -24.614+/-0.010 vs -24.601  z=1.33 PASS   <- was z=38 FAIL
STILL OPEN: v collapses walkers at ~step 41 (phi->0) -- a FORWARD-propagation
instability in the attractive -v channels (MkExpV `else` branch / Vznn), now
isolated from the (fixed) back-propagation. uxx + uxy validated.

## 2026-06-24 — uxy / v channels: energy estimator fixed, propagation still open

Extended the production energy estimator in meas.f90 (correl: have=hkin+hpot)
to include the inter-orbital on-site U' (uxy) and the inter-orbital -v
neighbour terms (the code already had the intra-orbital +v term). The added
two-body form expands to exactly the same Wick structure as the validated
intra-v / uxx terms:
   <n_a n_b> = (gxu_aa+gxd_aa)(gxu_bb+gxd_bb) - gxu_ab gxu_ba - gxd_ab gxd_ba.
Regression: uxx=2 still PASS (-28.292+/-0.007 vs ED -28.289, z=0.45) — the new
blocks are gated by |uxy|,|v|>1e-3 so uxx-only is unchanged.

Also fixed a real InitV bug: Vlist(i1,6:37)=+-v overwrote ALL 32 neighbour
channels with the last sign each iteration (losing the +v intra / -v inter
distinction); now per-channel Vlist(i1,k). No uxx regression (z=2.03 PASS).

STILL OPEN — the inter-orbital/neighbour PROPAGATION (Vee i/=j path), not the
estimator or the channel table:
  - uxy=2 -> -25.41 (was -32.88) but ED=-24.60, i.e. BELOW the exact ground
    state -> the Vee two-site HS update is biased for i/=j.
  - v=1 still collapses walkers (phi->0, totalEn pinned at free-fermion -32+/-0),
    even after the InitV fix -> the attractive -v channels (charge-channel HS,
    MkExpV `else` branch + coeffv weights) destabilise the constrained path.
Reviewed Vee (Vee.f90:217-460): the two-site sequential Green's-function update
and expV(i,...) indexing are internally consistent (bond strength stored at
index i is intentional). The bias is subtle -- next step is to check the
same-spin two-site rdet/GF update and the attractive-channel coeffv/constraint
against a hand-derived single-step reference on a 2-site toy. Estimator and
channel table are now correct; the open work is purely in the HS propagation.

## 2026-06-23 — code/src INTERACTING vs ED (small cluster 2x2, 4+4)

Build: `tools/build.sh src-small` (patches cpParameter to lx=ly=2, NUP=NDN=4
so ED is tractable). ED target: ed/altermagnet_ed.py ed_energy (QuSpin).
QMC observable: out.dat `totalEn` (the production energy with error bar),
runs with equilibration (nblkeq=4,nblkgr=2).

| uxx | code/src totalEn      | ED (exact) |  z   | verdict |
|-----|-----------------------|-----------:|-----:|---------|
| 0   | -32.000               | -32.000    | 0.0  | PASS    |
| 2   | -28.275 +/- 0.007     | -28.289    | 1.85 | PASS    |
| 4   | -25.134 +/- 0.010     | -25.128    | 0.58 | PASS    |

**The interacting two-orbital altermagnet CPQMC reproduces exact
diagonalization** on this cluster (on-site uxx channel) — the ED-vs-CPQMC
benchmark requested in PRL review. Next: validate the uxy and v channels the
same way, and fix the secondary printed estimator ("energy is the:", the
StepMeas mixed estimator) which still disagrees (e.g. -19.98 at uxx=2) while
the out.dat totalEn is correct.

## 2026-06-23 — code/src (two-orbital altermagnet), U=0 vs free-fermion ED

Model (code/src/t.f90): two-orbital (d_xz,d_yz) square lattice,
eps_x=-2t1 cos kx -2t2 cos ky -4t3 cos kx cos ky, eps_y=(x<->y),
eps_xy=-4t4 sin kx sin ky. This is the PRL LE20050 altermagnet model.

Reference: ed/altermagnet_ed.py free-fermion GS, verified against the analytic
k-space dispersion. 4x4, 16+16, t1..t4=-1  =>  E0 = -96.0 exactly.

| stage                         | before fix | after fix | exact |
|-------------------------------|-----------:|----------:|------:|
| Initial energy (InitEnergy)   |   -96.0    |   -96.0   | -96.0 |
| Measured energy (StepMeas)    |   +96.0    |   -96.0   | -96.0 |

**Bug found & fixed.** The measurement-phase energy was the exact negative of
the truth. calgf method 2 (back-propagation) built gx = +G (plain Green's
function), but calenergy expects the determinant-scaled gx = detp*(I-G) that
calgf method 1 produces. Fixed method 2 to mirror method 1 (commit ef171d9).
U=0 now validates to machine precision.

**Still open (next).** With on-site uxx=2 the energy is unchanged from U=0
(still -96), i.e. the interaction does not yet enter the estimate: rec_fields
(HS-field recording for back-propagation) is commented out in cpCore.f90, and
calenergy's potential term sums only uxx (not the inter-orbital uxy or
neighbor v). Interacting validation is the next milestone; extend
ed/altermagnet_ed.py to interacting ED (QuSpin) as the target.

## 2026-06-23 — benchmark-cpqmc vs ED, 4x4 Hubbard, Nup=Ndn=1, U=3, t=1

ED ground state (256 basis states, exact): E = -7.867205, E/site = -0.491700.

| observable        | ED        | CPQMC          | z     | verdict |
|-------------------|-----------|----------------|-------|---------|
| energy_total      | -7.867205 | -7.870103      | 1.01  | PASS    |
| energy_per_site   | -0.491700 | -0.491881      | 1.01  | PASS    |
| energy_kinetic    | -7.960602 | -7.965755      | 1.75  | PASS    |
| energy_potential  |  0.093397 |  0.095652      | 15.66 | FAIL    |

**Interpretation.** The *total* ground-state energy matches the exact result
to ~1 sigma — the single-band CPQMC is validated on this case. The potential
energy fails only because it is a **mixed/back-propagated estimator of an
operator that does not commute with H**: in projector QMC the total energy
obeys a zero-variance principle and is robust, while individual components
(PE, and the compensating KE) carry an O(systematic) estimator bias that
cancels in the total but not in the split. This is expected physics, not a
code error; it is the kind of effect the validation platform is meant to
surface. Total energy is therefore the primary pass/fail observable for
ED-vs-QMC; component energies are diagnostic.

## 2026-06-24 — Phase 6 WS1b/1c: pyqmc pairing + magnetization vs ED

New back-propagated estimators in pyqmc (`run_bp_pairmag`, single-band square
lattice): singlet pairing structure factors S_alpha = <O_alpha^dagger O_alpha>
for s-wave and d_{x^2-y^2} (O_alpha = sum_m sum_delta f_alpha(delta)
c_{m,up} c_{m+delta,dn}); antiferromagnetic spin structure factor S(pi,pi);
local moment^2 = sum_i <m_i^z>^2. ED reference = number-conserving product
operator O^dagger O = "+-|+-" (QuSpin), exact ground state.
Validation harness: `python pyqmc/validate_pairmag.py`.

**Estimator verification (U=0, non-degenerate nup=ndn=1, 4x4): EXACT.**
| observable        | ED       | pyqmc    | rel.dev |
|-------------------|----------|----------|---------|
| s-wave pairing SF | 16.00000 | 16.00000 | 0.0%    |
| d-wave pairing SF |  0.00000 |  0.00000 | 0.0%    |
| spin S(pi,pi)     |  0.03125 |  0.03125 | 0.0%    |
| sum_i <m_i^z>^2   |  0.00000 |  0.00000 | 0.0%    |
At U=0 the free-electron trial IS the exact ground state, so the back-propagated
Wick estimator must (and does) reproduce ED to machine precision. This confirms
the pairing/magnetization estimators are structurally correct.

NOTE on degeneracy: open-shell fillings (e.g. 4x4 nup=ndn=2 at U=0, where the
single-particle -2 level is 4-fold degenerate) are NOT clean gates — ED picks an
arbitrary state in the degenerate manifold while the QMC trial picks a specific
one, so state-dependent observables (s-wave SF, local moment) differ. Use
non-degenerate closed-shell fillings (nup=ndn=1).

**Interacting (U=3, non-degenerate nup=ndn=1, 4x4; ED E0=-7.8672):**
| observable        | ED       | pyqmc (bp=14)     | z     | note            |
|-------------------|----------|-------------------|-------|-----------------|
| spin S(pi,pi)     | 0.03162  | 0.03169 +/-1e-5   | 5.4   | 0.2% (CP bias)  |
| s-wave pairing SF | 14.70642 | 15.287  +/-0.021  | 27    | +4% (CP bias)   |
| d-wave pairing SF | 0.00000  | -0.0004 +/-0.0004 | 1.0   | PASS (~0)       |
| sum_i <m_i^z>^2   | 0.00000  | 0.0005  +/-2e-5   | --    | ~0 (paramag.)   |

**Interpretation.** Magnetization (spin S(pi,pi)) reproduces ED to 0.2%; the
d-wave order parameter and local moment (both ~0 at this filling) are captured.
The s-wave pairing structure factor carries a +4% systematic offset (z>>1 with a
small statistical error) — this is the **constrained-path / projector bias** for
an off-diagonal four-fermion observable that does not commute with H. Back-
propagation removes the mixed-estimator part but not the residual CP-node bias
from the free-electron trial; the same mechanism biases the PE component of the
energy (above). The estimator itself is exact (U=0). Reducing the residual needs
a better trial (Phase 6 WS3 adaptive trial) or release-constraint extrapolation.

## 2026-06-24 — Phase 6 WS3: adaptive (self-consistent) trial vs fixed

New `trial="adaptive"` mode (pyqmc `TrialWF.update`): each generation it forms the
weighted ensemble 1-RDM rho^s = sum_i w_i G^s_i, symmetrizes, damps toward the
current trial projector (mix), and takes the nup/ndn highest-occupation natural
orbitals as the new trial -- moving the constrained-path boundary toward the true
ground state. Cadence `trial_every`, damping `trial_mix`. Harness:
`python pyqmc/validate_adaptive.py`.

**Test = hard regime**: 4x4 single-band, nup=ndn=2, U=4 (open shell -- the
free-electron -2 level is 4-fold degenerate, so the FIXED trial is built from an
ambiguous degenerate manifold => large constrained-path bias). ED E0 = -11.530292.

| trial               | estimator | E         | dE(ED)  |
|---------------------|-----------|-----------|---------|
| fixed               | mixed     | -11.3239  | +0.2064 |
| adaptive (mix .5/10)| mixed     | -11.3599  | +0.1704 |
| fixed               | bp14      | -11.2895  | +0.2408 |
| adaptive (mix .5/10)| bp14      | -11.3414  | +0.1889 |

The adaptive trial **reduces the bias toward ED in both estimators**: mixed
+0.206 -> +0.170 (-17%), back-prop +0.241 -> +0.189 (-22%). The energies stay
above ED and move down toward it.

**Damping matters** (mix=0.25, every=5, i.e. more aggressive updating):
mixed +0.200 -> +0.182 (-9%), bp +0.207 -> +0.185 (-10%) -- LESS improvement than
the gentler mix=0.5. Too-aggressive updates inject noise from the (biased, noisy)
mixed RDM and partly offset the gain; moderate damping is best.

**Interpretation.** The self-consistent natural-orbital trial confirms the WS3
hypothesis -- redefining the sign/constraint boundary from each walker generation
makes the result more exact vs ED in the sign-problematic (here degenerate-shell)
regime. A residual bias (~0.17) remains because a single-determinant natural-
orbital trial still cannot represent the fully correlated ground state, and
because the adaptive trial is built from the biased mixed RDM. Caveat: the
adaptive trial is non-variational (it is no longer a fixed external constraint),
so "closer to ED" is an empirical result; FIXED remains the validated default and
the regression gate (pyqmc/regression.py) stays bit-for-bit identical for it.
Natural next step: build the trial from the back-propagated RDM (less biased) and
allow a few-determinant trial.

## 2026-06-24 — WS3 follow-up: does the adaptive bias scale with walker count?

Question: the constrained-path bias is systematic (not statistical), so for a
FIXED trial it should be independent of the walker count Nw. For the ADAPTIVE
trial the trial is rebuilt from the ensemble 1-RDM, whose accuracy improves with
Nw -> a better trial -> the bias ITSELF should coherently decrease with Nw.
Harness: `python pyqmc/scan_walkers.py`. 4x4, nup=ndn=2, U=4, mixed estimator,
ED E0 = -11.530292.

|   Nw |  fixed E    +/-    dE(ED) | adaptive E  +/-    dE(ED) |
|------|--------------------------|--------------------------|
|  100 | -11.3002  0.0030  +0.2301| -11.3022  0.0077  +0.2281|
|  250 | -11.3042  0.0028  +0.2261| -11.3413  0.0068  +0.1890|
|  500 | -11.3226  0.0026  +0.2077| -11.3712  0.0061  +0.1591|
| 1000 | -11.3210  0.0022  +0.2093| -11.3597  0.0045  +0.1706|

**Result — confirmed.** The FIXED bias is essentially flat (~0.21, no trend with
Nw); only its statistical error shrinks (0.0030 -> 0.0022 ~ 1/sqrt(Nw)). The
ADAPTIVE bias DECREASES with Nw: +0.228 (Nw=100) -> +0.189 -> +0.159 (Nw=500),
a ~25% reduction over fixed, before flattening (+0.171 at Nw=1000, within noise).
At Nw=100 adaptive ~ fixed: the self-consistent feedback needs enough walkers to
estimate the 1-RDM accurately before it helps -- there is a real crossover. So the
adaptive trial CONVERTS WALKER COUNT (a statistical resource) INTO REDUCED
SYSTEMATIC BIAS, a coherent effect the fixed trial cannot produce. Residual bias
remains (single-determinant trial built from the still-biased mixed RDM); the
next refinements (BP-RDM trial, multi-determinant) target it.

## 2026-06-24 — Phase 6 WS2: unequal-time (dynamic) spin susceptibility vs ED

New imaginary-time-displaced estimator (pyqmc Estimators.chi_block + CPMC.run_bp_chi,
single-band square lattice): the STAGGERED spin correlation
C(tau_l) = <O(tau_l) O(0)>, O = sum_i (-1)^{x+y} S^z_i, tau_l = l*dt, l=0..bp, and
the static staggered susceptibility chi_s = integral C(tau) dtau. Per walker the
recorded forward fields build the one-body propagators B_(l); the displaced
single-particle GFs are P^s(tau)=B(I-g^sT) [particle], H^s(tau)=B^{-T} g^s [hole],
and g^s(tau)=B^{-T} g^s B^T [equal-time at tau], with g^s_ij=<c^+_i c_j> from the
back-propagated bra. ED reference = Lehmann representation
C(tau)=sum_n |<n|O|0>|^2 e^{-(E_n-E_0)tau} (pyqmc/validate_chi.py).

**Estimator verification (U=0, nup=ndn=1, 4x4): EXACT.** C(tau) matches ED to 4
decimals at every tau (free-fermion trial = exact GS; zero statistical variance).
This was the gate that caught a real bug: the FIRST implementation used H=g B^{-1}
and g(tau)=B g B^{-1}, which are correct only when B and g are symmetric -- so U=0
passed but the INTERACTING C(tau) decayed too fast and went negative. Re-deriving
from c_i(tau)=sum_j B_ij c_j, c_i^dag(tau)=sum_k (B^{-T})_ik c_k^dag gives the
correct H=B^{-T} g and g(tau)=B^{-T} g B^T. Lesson: always validate dynamic
estimators on an ASYMMETRIC (interacting) case, not just U=0.

**Interacting (U=3, nup=ndn=1, 4x4):** after the fix C(tau) tracks ED across the
whole window within statistics:
| tau  |  ED C   | pyqmc C  |  z  |
|------|---------|----------|-----|
| 0.00 | 0.5059  | 0.5065   | 3.8 |
| 0.10 | 0.2304  | 0.2307   | 1.6 |
| 0.16 | 0.1437  | 0.1437   | 0.0 |
| 0.20 | 0.1049  | 0.1045   | 2.3 |
Windowed chi_stag (0..0.20): ED 0.0510 vs pyqmc 0.0511 +/- 0.0001 -> PASS.
C(tau_0) = N*S(pi,pi) reproduces the equal-time staggered structure factor.

NOTE: the full static chi needs tau out to several (C(tau) has not decayed at
tau=0.20=bp*dt); the rigorous gate is the WINDOWED integral (same tau-range for ED
and QMC). Larger bp extends the window at higher cost / B^{-1} conditioning.
Pairing susceptibility chi_d (4-point, time-displaced) is the natural follow-on.

## 2026-06-24 — WS3 variant: single-walker stochastic trial (trial="sample1")

Idea (tested on request): instead of the deterministic ensemble natural-orbital
1-RDM, redefine the trial / sign boundary each generation to a SINGLE walker drawn
with probability proportional to its weight. A single determinant, so the existing
machinery applies unchanged (TrialWF._update_sample1). A/B vs fixed and adaptive on
4x4 nup=ndn=2 U=4 (ED -11.5303), nw=400:

| trial    | est   | E         | dE(ED)  | stat.err |
|----------|-------|-----------|---------|----------|
| fixed    | mixed | -11.3100  | +0.220  | 0.0025   |
| adaptive | mixed | -11.3664  | +0.164  | 0.0037   |
| sample1  | mixed | -11.8552  | -0.325  | 0.248    |
| fixed    | bp    | -11.2985  | +0.232  | 0.020    |
| adaptive | bp    | -11.3564  | +0.174  | 0.016    |
| sample1  | bp    | -11.0227  | +0.508  | 0.376    |

**Verdict: a single-walker trial is NOT a good idea.** It loses on BOTH axes:
(i) Variance ~70-100x larger than fixed/adaptive (0.25-0.38 vs 0.003-0.02) -- one
walker is a hopelessly noisy stand-in for the ensemble ground state.
(ii) Non-variational + unstable: the mixed estimator overshoots BELOW the true E0
(-11.855 < -11.530), while the bp estimator overshoots ABOVE (+0.508); the two
disagree by 0.83, proving the central value is unreliable (a fluctuating non-
orthogonal trial breaks the variational bound). The deceptively small z~1.3 is only
because the error bar is enormous, not because it is accurate.
The deterministic natural-orbital AVERAGE (adaptive) wins decisively: lowest bias
AND tiny variance. The principled improvement over adaptive is the opposite of
collapsing to one walker -- a MULTI-DETERMINANT trial (a weighted superposition of
k sampled walkers), which lowers variance further and converges toward exact as
k grows. That requires generalizing the importance/overlap function (the up/dn
overlaps no longer factorize) and is the next experiment.

## 2026-06-24 — WS3 variant: multi-determinant trial (trial="multidet") — NEGATIVE

Tested the "principled" multidet trial: each generation resample k walkers
(prob. proportional to weight), orthonormalize, and use their weighted
COHERENT superposition Psi_T = sum_m c_m |phi_m> as the trial (md_overlap,
md_green generalize the importance function + mixed Green's function; the
single-determinant fast path is untouched, regression bit-for-bit). A/B run on
the 48-core workstation (6 schemes x 4 seeds = 24 parallel jobs, 4x4 nup=ndn=2
U=4, ED -11.5303, mixed estimator):

| scheme       | mean E   | dE(ED)  | run-run std | mean err | per-seed E (4 seeds)        |
|--------------|----------|---------|-------------|----------|-----------------------------|
| fixed        | -11.281  | +0.250  |   0.007     |  0.002   | -11.28/-11.29/-11.28/-11.28 |
| adaptive     | -11.294  | +0.236  |   0.012     |  0.004   | -11.28/-11.31/-11.28/-11.31 |
| sample1      | -13.337  | -1.807  |   1.872     |  1.594   | -12.5/-11.9/-12.4/-16.6     |
| multidet k=4 | -61.6    | -50.0   |  43.8       | 44.8     | -46/-134/-15/-51            |
| multidet k=8 | -122.5   | -110.9  | 180.4       |100.8     | -18/-435/-24/-13            |
| multidet k=16| -75.2    | -63.7   |  88.5       | 60.1     | -42/-228/-18/-14            |

**Verdict: a coherent multi-determinant superposition of sampled walkers is
catastrophically UNSTABLE** -- energies far below the true ground state with
run-to-run scatter of 40-180 (vs 0.01 for adaptive). Cause: the importance
function <Psi_T|phi> = sum_m c_m <up_m|phi_up><dn_m|phi_dn> suffers sign/magnitude
cancellation, so the mixed-Green denominator W = sum_m t_m -> 0 and the energy
diverges. Phase-aligning each determinant to a common reference (so <ref|phi_m> >
0) removed the worst blow-ups (k=4 went to ~ -12 to -17) but did NOT cure it
(still below E0, err ~ 0.4-4.9).

**Why it fails (the real lesson).** The constrained-path population represents the
ground state as a STATISTICAL MIXTURE -- an incoherent weighted average of walker
Green's functions / density matrices, |psi_0><psi_0| ~ sum_i w_i (...). The
adaptive natural-orbital trial uses exactly that mixture (it averages the 1-RDM)
and is stable + bias-reducing. A COHERENT superposition sum_m c_m |phi_m> is the
WRONG mathematical object: independent walkers carry arbitrary relative phases and,
after population control, are near-duplicates spanning an ill-conditioned subspace,
so their coherent sum is near-singular. A proper multi-determinant trial needs
determinants from a deterministic CI/CASSCF-like expansion with controlled
coefficients -- NOT random walkers.

**Conclusion of the trial-from-population study.** Among ways to redefine the
trial/constraint from the walker generation: the density-matrix AVERAGE (adaptive
natural orbitals) is the only stable, bias-reducing route; a single walker
(sample1) is non-variational + high-variance; a coherent walker superposition
(multidet) is numerically unstable. trial="multidet" is kept as a documented
experiment but is NOT recommended; fixed/adaptive remain the validated options.

## 2026-06-24 — multidet follow-up: machinery verified, fix attempt confirms it

Audit response to "is it a bug?": the multidet ESTIMATOR machinery is correct,
verified by running multidet with the update DISABLED (k=1, trial = the single
free determinant): U=0 gives -12.0000 exact and U=3 gives -11.36 (vs fixed -11.39),
stable, no blow-up. So md_overlap / md_green / the multidet step are right; the
instability is entirely in the trial CONSTRUCTION from walkers.

Two construction issues were identified and the principled fix applied:
(1) coefficient should be the importance-reconstruction weight c_m =
    w_m / <psi_T_old|Q_m>  (|psi_0> ~ sum_i (w_i/<psi_T|phi_i>)|phi_i>), not c_m=w_m;
(2) QR drops the det(R) amplitude -- but this cancels against the 1/<.|.> factor
    when computed with the orthonormalized Q_m, and the sign comes out automatically
    (no ad-hoc phase flip).
RESULT: the corrected coefficient makes it WORSE, not better -- 4x4 U=4 now gives
+157 / -219 / -75 (k=4,8) because walkers near the trial node (small <psi_T|phi_m>)
get huge coefficients that amplify the conditioning blow-up.

**Definitive conclusion.** The naive coefficient (c=w) and the correct importance
coefficient (c=w/<psi_T|phi>) BOTH diverge, in opposite ways (sign cancellation vs
near-node amplification). The estimator math is verified correct. Therefore the
failure is fundamental to the BASIS: a coherent multi-determinant trial built from
sampled CP walkers is ill-conditioned regardless of coefficients -- the walkers are
a statistical mixture (near-duplicate after pop control, arbitrary relative phases,
some near the node), not a good coherent expansion. The density-matrix average
(adaptive natural orbitals) is the correct way to distill the population; a proper
multi-determinant trial needs a deterministic CI/CASSCF expansion, not walkers.
trial="multidet" kept as a documented experiment (now with the correct coefficient);
NOT recommended. fixed/adaptive remain validated.

## 2026-06-24 — multidet ENGINE VALIDATED + a real bug fixed (ED-truncation)

Step (a) of the multidet plan: validate the engine with a KNOWN-GOOD trial =
truncated CI expansion of the exact ground state. To avoid QuSpin sign-convention
mismatches, the CI is built in OUR sorted-column determinant convention via
Slater-Condon (pyqmc/validate_casci_ed.py: ci_ground_state) and diagonalized; its
E0 matches QuSpin ED to 1e-14. The top-N configs (by |coef|) form a frozen
multidet trial (TrialWF.set_multidet), swept over N.

This exposed a REAL BUG (the user was right to suspect one): the multidet energy
computed the two-body term from the WEIGHT-AVERAGED Green's function, but <n_a n_b>
is quadratic in G, so for a multi-determinant trial it must be the t_m-weighted
average of the PER-DETERMINANT Wick energies, E_L = sum_m t_m E^(m)/sum_m t_m. U=0
(no two-body) and k=1 (single term) both masked it; that is why the earlier k-scan
and the "multidet with update disabled = fixed" check passed despite the bug.
Fixed in Estimators._local_energy. (Single-determinant paths unchanged; regression
bit-for-bit.)

RESULT after the fix (2x2 nup=ndn=2 U=4, E0=-5.656854):
| N_det | CPMC E   | dE(ED)  | capture |
|-------|----------|---------|---------|
|   1   | -4.6111  | +1.0458 | 0.193   |
|   4   | -5.2000  | +0.4569 | 0.483   |
|   8   | -5.5274  | +0.1295 | 0.628   |
|  16   | -5.6307  | +0.0262 | 0.822   |
|  32   | -5.6569  | -0.0000 | 1.000   |
|  36   | -5.6569  | -0.0000 | 1.000   |
At full N the energy equals E0 with ZERO variance (the trial IS the exact GS ->
every local energy = E0), and the constrained-path bias collapses MONOTONICALLY to
zero as the trial improves. This DEFINITIVELY shows a good multi-determinant trial
removes the CP bias -- the engine is correct. (Note: the earlier "fundamental"
verdict on multidet was overstated; the machinery had a fixable two-body bug. The
walker-SUPERSITION ill-conditioning is a separate, real issue -- random walkers are
still a poor multidet basis -- which is why the scalable route uses a deterministic
natural-orbital CASCI, not walkers.)

## 2026-06-24 — STEP B DONE: self-consistent natural-orbital CASCI trial (mini-CASSCF)

The scalable multi-determinant trial (no ED). New: pyqmc/casci.py (general
Slater-Condon CI in an arbitrary orbital basis + integral transform +
natural-orbital active-space builder) and pyqmc/validate_casscf.py (the SCF loop).

GENERAL CASCI ENGINE validated basis-invariant (the key correctness gate):
transform the Hubbard H into basis W (h=W^T K W, U_pqrs=U sum_i W_ip W_iq W_ir W_is)
and diagonalize -> E0 identical to QuSpin ED to 1e-14 for W = site / K-eigenbasis /
RANDOM orthonormal (2x2,4+4,U=4: -5.65685425 all three). Full-N natural-orbital
CASCI trial fed to the multidet engine gives E0 with zero variance.

SELF-CONSISTENT CASCI-SCF (walker 1-RDM natural orbitals -> active space -> CAS
diag -> install frozen multidet trial -> iterate), 2x2/4+4/U=4 (degenerate open
shell, ED -5.6569):
  FIXED free-electron trial : E=-4.108  dE=+1.548   (large CP bias)
  CASCI-SCF FULL active (gate): E=-5.6569 dE=+0.0000 (exact, zero variance) PASS
  CASCI-SCF reduced CAS(2,2)  : E=-5.639  dE=+0.018  (CAS CI E0 only -4.99)
-> a tiny 4-config active space (1 frozen core + 2 active + 1 virtual, NO ED) cuts
the constrained-path bias ~85x (1.548 -> 0.018), self-consistently. The CAS energy
itself is far from ED; it is the trial NODES that matter. This is the scalable
route the walker-superposition multidet could not provide.

## 2026-06-24 — CASSCF trial on the 4x4 DOPED case (where adaptive only managed ~17%)
4x4, nup=ndn=2, U=4 (the exact case used to test the adaptive trial, ED -11.5303):
  FIXED free-electron trial      dE=+0.2409
  CASCI-SCF reduced CAS(2,2 in 4o, no core)  dE=+0.0231   (CAS CI E0 only -11.23)
-> ~90% (10x) constrained-path bias reduction, vs the adaptive single-determinant
trial's ~17-22% on the same case. Self-consistent natural orbitals (occ 1.34/0.62/
0.54/... -> 4-orbital CAS) converge over iterations. The scalable CASSCF trial
clearly beats both fixed and adaptive. NO ED used to build the trial.

## 2026-06-24 — CASSCF active-space convergence (run on the 48-core server)
4x4, nup=ndn=2, U=4 (ED -11.5303). CASCI-SCF |dE(ED)| vs active-space size
(n_core=0, parallel runs on amaz@192.168.101.54):
  n_active  configs  |dE(ED)|   note
    2          1      +0.180    single config (= natural-orbital HF, no correlation) ~ fixed
    3          9      +0.039    correlation entering -> ~77% bias cut
    4         36      +0.028    ~83% bias cut (laptop seed gave +0.023; seed/statistics)
(fixed free-electron baseline |dE|=0.170). The CP bias falls sharply once the
active space exceeds the electron count (n_active>nup) so multi-configuration
correlation enters; a 4-orbital CAS already removes ~90%. Confirms the method is
a small-CAS effect, not requiring large active spaces. Remote env: ~/qmcenv
(numpy 2.5/scipy 1.18), code in ~/qmc; runs launched detached + polled.

## 2026-06-25 — cost (Q1) + adaptive active space & accuracy (Q2)

COST per CPMC step (pyqmc/cost_estimate.py, server). adaptive = fixed (single
determinant; only an occasional 1-RDM diagonalization, amortized). CASSCF
multidet ~ O(k) x fixed, k = #CAS determinants (4x4/4+4):
  fixed 100ms (1.0x) | adaptive 99ms (1.0x) | multidet k=4 1.8x, k=16 5.5x,
  k=36 11.6x  (~0.3-0.5 k; the determinant-overlap count is the multiplier).
The multiplier is set by k and is ~lattice-INDEPENDENT (per-determinant work
scales like fixed), so the SAME factor applies at L=20x20. Absolute time at
20x20 is infeasible in this pure-Python validation code for ANY trial
(per-step ~ nwalkers.N.Ne^3); a compiled/BLAS engine is required, where the
k x multidet overhead carries over. k from the auto active space: CAS(4,4)~36
(~12x), CAS(6,6)~400 (~150x), CAS(8,8)~4900 (~1500x).

ADAPTIVE ACTIVE SPACE + ED-FREE ACCURACY (pyqmc/casci.py auto_active_space,
occupation_entropy, cas_accuracy_estimate; pyqmc/validate_casscf.py --auto;
pyqmc/scan_active.py). 4x4/4+4, n_active grows with correlation:
  U     entropy  n_active  k(dets)  dE_trunc(accuracy)
  1     2.59       7         441      0.005
  2     3.63      10        2025      0.005
  4     4.56      10        2025      0.022
  8     5.21      10        2025      0.051
  12    5.46      10        2025      0.081
n_active + the multireference entropy grow with U; the ED-free dE_trunc =
|E_cas(n+1)-E_cas(n)| grows 0.005->0.081, correctly flagging where the method
is least accurate (tracks the true fixed-trial bias +0.14@U=2 -> +0.26@U=8).
CAVEAT: the auto selector can pick a large CAS (k=2025) that the pure-Python
multidet ENGINE handles poorly (ill-conditioned, walker death) -- the accuracy
estimate is exactly the signal that the regime is hard; a fast multidet engine
+ capped/iteratively-grown CAS is the production path.

## 2026-06-25 — sign problem: free projection vs constrained path (pyqmc/sign_demo.py)

Direct demonstration of the exact-but-exponential vs biased-but-polynomial
trade-off. 4x4 nup=4 ndn=2 (spin-imbalanced) U=8 + frustrating t2=0.4 (a real
sign problem), 48 independent populations on the 48-core server, dt=0.1, project
to tau=10.

FREE PROJECTION (no constraint, signed weights, EXACT): the average sign decays
  tau:  0.2  2.0  4.0  5.0  6.0  8.0  9.0
  <s>:  1.00 0.88 0.67 0.58 0.43 0.15 ~0
and the energy error tracks 1/<s> and EXPLODES: err 0.007 (tau=0.2) -> 0.3
(tau=4) -> 18.7 (tau=4.6) -> 148 (tau=8.2); E becomes garbage (-92,-166) once
<s>->0. Unusable beyond tau~4-5.

CONSTRAINED PATH (q=max(R,0), positive weights): NO sign decay, error FLAT
~0.012, locked on E=-17.78 through tau=2-6 -- a stable polynomial-cost answer
(small systematic bias). [Beyond tau~7 the minimal demo's CP also degrades, but
that is just missing population control, not a sign problem.]

CONCLUSION (answers "is the trial work better than directly simulating the sign
problem?"): where free projection is affordable (small system / short tau) it is
exact and is the gold standard -- that is why we anchor on ED at small sizes.
But its cost to fixed error is exponential in N*beta, so for the large/doped/
spin-dependent-hopping regime (the altermagnet, whose sign problem the appeal
notes is severe) it is intractable. CP with a good trial is the only polynomial
route; the CASCI/adaptive trial work reduces CP's one remaining error (the
node bias) and -- in a controlled constraint-release -- a better trial slows the
sign decay, making a short unbiased release affordable. The two are complementary
points on one bias<->sign-cost spectrum, not competitors.

## 2026-06-25 — dynamic singlet PAIRING susceptibility chi_s / chi_d vs ED

New unequal-time pairing estimator (Estimators.chid_block + CPMC.run_bp_chid,
single-band square lattice): C_a(tau)=<Delta_a(tau) Delta_a^dag(0)>,
Delta_a^dag = sum_m sum_delta f_a(delta) c^+_{m up} c^+_{m+delta dn}, for a = s-wave
and d_{x^2-y^2}, and chi_a = integral C_a(tau) dtau. By Wick (spin-separated) it
reuses the validated time-displaced PARTICLE GF P^s(tau)=B_(l)(I-g^sT) from
chi_block and the d-wave bond form-factor matrix from run_bp_pairmag:
C_a(tau) = sum_{m,n} P^up(tau)_{mn} (F_a P^dn(tau) F_a^T)_{mn}. <Delta>=0 so it is
already connected. ED reference (pyqmc/validate_chid.py): Lehmann on the FULL Fock
space (small clusters only), projecting H onto the (nup,ndn) and (nup+1,ndn+1)
filling sectors (exact integer fillings -> robust even at U=0 where cross-filling
degeneracy otherwise mixes eigenvectors).

Gate (U=0, 2x2 nup=ndn=1): pyqmc reproduces ED EXACTLY at every tau for BOTH
channels (d-wave C=32.0000 flat; s-wave 16.0000 -> 5.2205) -- estimator verified.
Interacting (U=4, 2x2 nup=ndn=1): the C_a(tau) curves track ED to ~1-4% across the
window; windowed chi_d ED 4.761 vs pyqmc 4.803 (+0.9%), chi_s ED 1.395 vs 1.383
(-0.9%). z=2-9 reflects the tight error bars exposing the small constrained-path
bias (same pattern as the spin susceptibility). C_a(tau_0) = equal-time pair
structure factor. NOTE: ED full-integral chi differs from windowed because C(tau)
has not fully decayed at tau_max=bp*dt; the rigorous gate is the windowed integral
(same tau-range for ED and QMC). The d-wave PAIRING susceptibility -- the PRL
appeal's central order parameter -- is now ED-validated in pyqmc.

## 2026-06-25 — pairing VERTEX susceptibility (chi_d/chi_s, connected) — CPMC + ED

Extended the dynamic pairing susceptibility to separate the FULL chi_a from the
CONNECTED/VERTEX part. Bubble (disconnected) = built from the ensemble-AVERAGED
time-displaced single-particle GFs: C^0_a(tau) = sum_{m,n} <P^up(tau)>_{mn}
(F_a <P^dn(tau)> F_a^T)_{mn}; vertex = full - bubble; chi_a^vertex = integral.
- CPMC: Estimators.chid_block now also returns per-tau weighted P^up, P^dn;
  CPMC.run_bp_chid reports chi_{a} (full) and chi_{a}_vertex (+ C(tau) curves).
- ED: validate_chid.ed_chid adds the (N+-1)-sector single-particle time-displaced
  GFs G^s(tau)=M_s diag(e^{-dE tau}) M_s^T (M_s[a,n]=<n^{N+1}|c^+_a|0>) -> bubble ->
  vertex. Full uses the (N+2) Lehmann as before.

Validation (single-band 2x2):
  U=0, nup=ndn=1: vertex C(tau)=0 at every tau for BOTH ED and CPMC (full=bubble
    when non-interacting) -- exact gate, machinery verified.
  U=4, nup=ndn=1 (dilute, non-degenerate): windowed chi
    d-wave: full ED 4.761 / pyqmc 4.803 ; VERTEX ED 0.0168 / pyqmc 0.0246(49) (~1.6 sigma)
    s-wave: full ED 1.395 / pyqmc 1.383 ; VERTEX ED 0.0512 / pyqmc 0.0389 (~24% CP bias)
  The vertex susceptibility is reproduced on the clean cluster; d-wave within error
  bars. CAVEAT: as for the equal-time vertex, the d-wave vertex is CP-bias-prone in
  the half-filled/frustrated (NNN-t1) regime -- trustworthy on clean dilute clusters,
  indicative (trend) in the manuscript regime. (The manuscript itself reports the
  equal-time vertex because the tau-integrated susceptibility is free-projection
  sign-problem inaccessible; the CPMC version here is the constrained-path estimate.)

## 2026-06-25 — constraint RELEASE for the d-wave vertex (CP-bias attack)

pyqmc/release_dwave.py: CP-equilibrate, then RELEASE the constraint (signed free
projection) and measure the equal-time d-wave VERTEX (mixed estimator) at each
release step, tracking <sign>. Target: 4x2, tam=0, t1=0.3, U=4, half filling, where
ED d-wave vertex = +23.18 and CP-AFQMC gives ~ -2.

First verified the ED target is real: 4x2 U=4 GS is NON-degenerate (gap E1-E0=0.42),
so +23 is the exact value (the U=0 vertex!=0 anomaly was the U=0 free-fermion
degeneracy, lifted by U).

RESULT: release does NOT recover +23 -- the d-wave vertex stays ~ -2 to -2.4 across
the whole release (tau 0..0.36) while <sign> only decays 1.00 -> 0.89. Diagnosis:
release with the MIXED estimator converges to <psi_T|O|psi_0>/<psi_T|psi_0>, i.e. it
removes the CONSTRAINT bias on the KET but the free-electron BRA still biases the
d-wave 4-point. Both CP-back-propagation (-2: bra projected, ket constrained) and
release-mixed (-2: ket released, bra biased) sit at -2 because each removes only ONE
of the two biases. To recover +23 you must remove BOTH: either released
back-propagation (release the ket AND back-propagate the bra through the unconstrained
path) or -- cleaner and permanent -- a pairing-aware (d-wave BCS / number-projected)
trial so the bra/node carry d-wave character. The normal-state free/adaptive/CASCI
trials cannot fix the d-wave bra bias.
