#!/usr/bin/env python3
"""Constrained-path AFQMC (CPMC) — Python port (single-band Hubbard + two-orbital
altermagnet).

Phase 5 of the validation platform: a from-scratch Python CPMC whose results are
checked against the ED reference (ed/hubbard_ed.py, ed/altermagnet_ed.py) and the
Fortran benchmark. Validated on the single-band case (4x4, 1up/1dn, U=3, ED
= -7.8672) and the two-orbital altermagnet on all three interaction channels
(uxx/uxy/v) plus equal-time correlations.

Algorithm (Zhang-Krakauer CPMC, ground-state projection):
  H = -t sum_<ij>,s c^+_is c_js + (density-density interactions: uxx/uxy/v)
  Trotter:  e^{-dt H} ~ e^{-dt K/2} e^{-dt V} e^{-dt K/2}
  Discrete Hirsch HS for each density-density term; importance sampling of the
  HS fields with the trial overlap ratio; constrained path (drop walkers whose
  overlap with the trial turns non-positive). Mixed + back-propagated estimators.

Phase 6 OOP structure (numerically identical to the pre-refactor monolith; gated
bit-for-bit by pyqmc/regression.py):
  LatticeModel   - K, e^{-dt K/2}, interaction term list + HS factors (problem def)
  TrialWF        - trial Slater dets defining the constraint/importance function;
                   `update()` hook is the seam for the adaptive trial (WS3)
  WalkerEnsemble - walker dets, weights, overlaps; reorthogonalize + pop control
  Propagator     - one CPMC step (importance-sampled HS), with optional recording
  Estimators     - mixed energy, back-propagated energy + correlations
  CPMC           - orchestrator wiring the above; run / run_bp / run_bp_obs
"""
from __future__ import annotations
import numpy as np


def square_hopping(lx, ly, t=1.0):
    n = lx * ly
    K = np.zeros((n, n))
    def idx(x, y): return (x % lx) * ly + (y % ly)
    for x in range(lx):
        for y in range(ly):
            i = idx(x, y)
            for (jx, jy) in (((x + 1) % lx, y), (x, (y + 1) % ly)):
                j = idx(jx, jy)
                K[i, j] -= t
                K[j, i] -= t
    return K


def am_hopping(lx, ly, t0=1.0, tam=0.0, t1=0.0, tp=0.0):
    """Single-band ALTERMAGNET spin-dependent hopping (the manuscript model,
    mc2duph.f90). Two anisotropy channels:
      tam (NN): up hops t0-tam on x, t0+tam on y; down rotated 90 deg. This gives a
        d_{x^2-y^2} spin splitting eps_up-eps_dn ~ (cos kx - cos ky) BUT PRESERVES
        (pi,pi) nesting -> does not, on its own, suppress AFM.
      t1 (NNN, diagonal): spin-dependent d_xy hopping, up = +t1 on (+-1,+-1) main
        diagonal and -t1 on the anti-diagonal; down opposite. This adds cos(kx+-ky)
        terms that BREAK (pi,pi) nesting -> the nesting-disrupting anisotropy the
        manuscript invokes.
    tam=t1=0 -> isotropic (= square_hopping(t0)). Convention K = -hopping (H=sum K c+c).
    Returns (K_up, K_dn)."""
    n = lx * ly
    Ku = np.zeros((n, n)); Kd = np.zeros((n, n))
    def idx(x, y): return (x % lx) * ly + (y % ly)
    for x in range(lx):
        for y in range(ly):
            i = idx(x, y); jx = idx(x + 1, y); jy = idx(x, y + 1)
            Ku[i, jx] -= (t0 - tam); Ku[jx, i] -= (t0 - tam)   # up: weak x
            Ku[i, jy] -= (t0 + tam); Ku[jy, i] -= (t0 + tam)   # up: strong y
            Kd[i, jx] -= (t0 + tam); Kd[jx, i] -= (t0 + tam)   # dn: strong x
            Kd[i, jy] -= (t0 - tam); Kd[jy, i] -= (t0 - tam)   # dn: weak y
            if t1 != 0.0:
                md1 = idx(x + 1, y + 1); md2 = idx(x - 1, y - 1)   # main diagonal
                ad1 = idx(x - 1, y + 1); ad2 = idx(x + 1, y - 1)   # anti-diagonal
                for j in (md1, md2):       # up:-t1 main, dn:+t1 (K = -tk_Fortran)
                    Ku[i, j] -= t1; Kd[i, j] += t1
                for j in (ad1, ad2):       # up:+t1 anti, dn:-t1
                    Ku[i, j] += t1; Kd[i, j] -= t1
            if tp != 0.0:                  # ISOTROPIC (spin-independent) NNN t' --
                md1 = idx(x + 1, y + 1); md2 = idx(x - 1, y - 1)   # standard t' on
                ad1 = idx(x - 1, y + 1); ad2 = idx(x + 1, y - 1)   # ALL diagonals,
                for j in (md1, md2, ad1, ad2):                     # SAME for both spins
                    Ku[i, j] -= tp; Kd[i, j] -= tp                 # (no spin/d-wave structure)
    return Ku, Kd


# ---- generic Slater-determinant overlaps / Green's functions ----
def _ov_spin(left, phi):
    """<left|phi> single-spin overlap determinant det(left^T phi)."""
    return np.linalg.det(left.T @ phi)


def _green(left, phi):
    """G_ij = <c^+_j c_i> = [phi (left^T phi)^{-1} left^T]_ij for bra `left`."""
    A = left.T @ phi
    return phi @ np.linalg.solve(A, left.T)


class LatticeModel:
    """Problem definition: hopping K, half-step propagator e^{-dt K/2}, and the
    interaction as a general density-density term list V n_{a,sa} n_{b,sb} with the
    discrete-HS factors (covers on-site uxx, inter-orbital uxy, neighbour v)."""

    def __init__(self, K, lx, ly, dt, U, uxy, v, K_dn=None):
        # SPIN-DEPENDENT hopping: K is the up-spin one-body matrix; K_dn (optional)
        # the down-spin one. K_dn=None => spin-independent (K_dn=K), so the existing
        # single-K behaviour is reproduced bit-for-bit. self.K / self.expK stay the
        # UP-spin objects (back-compat); _dn variants drive the down channel.
        self.K = np.asarray(K, float)
        self.K_dn = self.K if K_dn is None else np.asarray(K_dn, float)
        self.n = self.K.shape[0]
        self.lx, self.ly, self.dt = lx, ly, dt
        self.U, self.uxy, self.v = U, uxy, v
        wu, vu = np.linalg.eigh(self.K)          # note: not 'v' (interaction param)
        wd, vd = np.linalg.eigh(self.K_dn)
        self.eigvecs = vu; self.eigvecs_dn = vd
        self.expK = vu @ np.diag(np.exp(-0.5 * dt * wu)) @ vu.T        # up e^{-dt K/2}
        self.expK_inv = vu @ np.diag(np.exp(+0.5 * dt * wu)) @ vu.T
        self.expK_dn = vd @ np.diag(np.exp(-0.5 * dt * wd)) @ vd.T     # dn e^{-dt K/2}
        self.expK_dn_inv = vd @ np.diag(np.exp(+0.5 * dt * wd)) @ vd.T
        self.terms = []       # (a, sa, b, sb, e1, e2, cf) for propagation
        self.pot_terms = []   # (a, sa, b, sb, V) for the energy estimator
        self._build_terms(U, uxy, v)

    def _hs(self, V):
        """Discrete-HS factors for exp(-dt V n_a n_b), field s in {+1,-1} (idx 0,1).
        Repulsive (V>0): spin channel, op1 e^{+a s - dtV/2}, op2 e^{-a s - dtV/2}.
        Attractive (V<0): charge channel, both e^{+a s - dtV/2}, weight e^{-a s + dtV/2}.
        cosh(a) = e^{dt|V|/2}."""
        dt = self.dt
        a = np.arccosh(np.exp(0.5 * dt * abs(V)))
        if V > 0:
            e1 = np.array([np.exp(a - 0.5 * dt * V), np.exp(-a - 0.5 * dt * V)])
            e2 = np.array([np.exp(-a - 0.5 * dt * V), np.exp(a - 0.5 * dt * V)])
            cf = np.array([1.0, 1.0])
        else:
            e1 = np.array([np.exp(a - 0.5 * dt * V), np.exp(-a - 0.5 * dt * V)])
            e2 = e1.copy()
            cf = np.array([np.exp(-a + 0.5 * dt * V), np.exp(a + 0.5 * dt * V)])
        return e1, e2, cf

    def _build_terms(self, U, uxy, v):
        # on-site Hubbard U on every (orbital-)site: n_{i,up} n_{i,dn}
        if U != 0:
            e1, e2, cf = self._hs(U)
            for i in range(self.n):
                self.terms.append((i, 0, i, 1, e1, e2, cf))
                self.pot_terms.append((i, 0, i, 1, U))
        if uxy == 0 and v == 0:
            return
        lxy = self.n // 2
        lx, ly = self.lx, self.ly
        def idx0(ix, iy, orb):
            return orb * lxy + (ix % lx) * ly + (iy % ly)
        # inter-orbital U' (uxy): n_r n_{r+lxy}, all 4 spin combos
        if uxy != 0:
            e1, e2, cf = self._hs(uxy)
            for r in range(lxy):
                for sa in (0, 1):
                    for sb in (0, 1):
                        self.terms.append((r, sa, r + lxy, sb, e1, e2, cf))
                        self.pot_terms.append((r, sa, r + lxy, sb, uxy))
        # neighbour V (+v intra-orbital, -v inter-orbital), +x and +y bonds
        if v != 0:
            hp = self._hs(v); hm = self._hs(-v)
            for ix in range(lx):
                for iy in range(ly):
                    for (jx, jy) in (((ix + 1) % lx, iy), (ix, (iy + 1) % ly)):
                        if (jx, jy) == (ix, iy):
                            continue                       # skip self-bond (e.g. ly=1)
                        for oi in (0, 1):
                            for oj in (0, 1):
                                a = idx0(ix, iy, oi); b = idx0(jx, jy, oj)
                                Vc = v if oi == oj else -v
                                e1, e2, cf = hp if oi == oj else hm
                                for sa in (0, 1):
                                    for sb in (0, 1):
                                        self.terms.append((a, sa, b, sb, e1, e2, cf))
                                        self.pot_terms.append((a, sa, b, sb, Vc))


class TrialWF:
    """Trial wavefunction defining the constrained-path boundary + importance
    function.

    mode="fixed" (default): the free-electron ground state (lowest nup/ndn
    orbitals of K), held constant -- the standard constrained-path trial.

    mode="adaptive" (Phase 6 WS3): a self-consistent trial. `update()` rebuilds
    the trial each generation from the current walker ensemble: it forms the
    weighted ensemble 1-body reduced density matrix rho^s = sum_i w_i G^s_i,
    symmetrizes it, damps it toward the current trial projector (mix), and takes
    the nup/ndn highest-occupation natural orbitals as the new trial. This moves
    the constraint boundary toward the true ground state, which should reduce the
    constrained-path / mixed-estimator bias where a fixed free-electron trial is
    poor (degenerate shells, doping / sign-problematic regime)."""

    def __init__(self, up, dn, mode="fixed", mix=0.5, rng=None, k=8):
        self.up = up
        self.dn = dn
        self.mode = mode
        self.mix = mix
        self.rng = rng
        self.k = k                                  # multidet: target #determinants
        self.nup = up.shape[1]
        self.ndn = dn.shape[1]
        self.multi = (mode == "multidet")
        self.frozen = False                         # if True, update() is a no-op
        # multi-determinant trial: list of determinant pairs + coefficients.
        # Starts as the single free-electron determinant (k=1); update() resamples.
        self.ups = [up]; self.dns = [dn]; self.coef = np.array([1.0])

    def overlap(self, phi_up, phi_dn):
        if self.multi:
            return self.md_overlap(phi_up, phi_dn)
        return _ov_spin(self.up, phi_up) * _ov_spin(self.dn, phi_dn)

    # ---- multi-determinant trial machinery ----
    def md_overlap(self, pu, pd):
        s = 0.0
        for m in range(len(self.coef)):
            s += self.coef[m] * _ov_spin(self.ups[m], pu) * _ov_spin(self.dns[m], pd)
        return s

    def md_green(self, pu, pd):
        """Multi-determinant mixed Green's functions G^s = sum_m t_m G^(m)_s / sum t_m,
        t_m = c_m <up_m|pu><dn_m|pd>. Returns (Gu, Gd, total_overlap)."""
        n = pu.shape[0]
        Gu = np.zeros((n, n)); Gd = np.zeros((n, n)); W = 0.0
        for m in range(len(self.coef)):
            ou = _ov_spin(self.ups[m], pu); od = _ov_spin(self.dns[m], pd)
            t = self.coef[m] * ou * od
            if t == 0.0:
                continue
            Gu += t * _green(self.ups[m], pu); Gd += t * _green(self.dns[m], pd); W += t
        return (Gu / W, Gd / W, W) if W != 0 else (Gu, Gd, 0.0)

    def _update_multidet(self, walkers):
        """Resample the trial as a k-determinant superposition of walkers drawn with
        probability proportional to weight. Importance-sampled CPMC represents the
        ground state as |psi_0> ~ sum_i (w_i / <psi_T|phi_i>) |phi_i>, so the COHERENT
        trial coefficient is c_m = w_m / <psi_T_old|Q_m> (Q_m = orthonormalized
        walker). This (a) restores the missing 1/<psi_T|phi> importance factor and
        (b) carries the correct relative sign automatically -- the det(R) amplitude
        dropped by QR cancels against the same factor in the overlap, so no ad-hoc
        phase alignment is needed. Coefficients normalized to O(1)."""
        w = np.clip(walkers.w, 0.0, None)
        if w.sum() <= 0:
            return
        kk = min(self.k, int((w > 0).sum()))
        idx = self.rng.choice(walkers.nw, size=kk, replace=False, p=w / w.sum())
        ups, dns, coef = [], [], []
        for j in idx:
            qu, _ = np.linalg.qr(walkers.phi_up[j])
            qd, _ = np.linalg.qr(walkers.phi_dn[j])
            olp = self.overlap(qu, qd)            # <psi_T_old | Q_m>  (old trial)
            if olp == 0.0:
                continue
            ups.append(qu); dns.append(qd); coef.append(w[j] / olp)
        if not ups:
            return
        self.ups = ups; self.dns = dns
        c = np.array(coef)
        self.coef = c / np.abs(c).max()           # normalize to O(1)
        self._refresh_olp(walkers)

    def _refresh_olp(self, walkers):
        walkers.olp = np.array([self.overlap(walkers.phi_up[i], walkers.phi_dn[i])
                                for i in range(walkers.nw)])

    def set_multidet(self, ups, dns, coef):
        """Install a fixed, externally-built multi-determinant trial (e.g. a
        CASCI / truncated-CI expansion) and freeze it (no resampling)."""
        self.multi = True
        self.mode = "multidet"
        self.frozen = True
        self.ups = [np.asarray(u, float) for u in ups]
        self.dns = [np.asarray(d, float) for d in dns]
        self.coef = np.asarray(coef, float)

    def update(self, walkers):
        if self.frozen or self.mode == "fixed":
            return
        if self.mode == "sample1":
            return self._update_sample1(walkers)
        if self.mode == "multidet":
            return self._update_multidet(walkers)
        if self.mode != "adaptive":
            return
        n = self.up.shape[0]
        ru = np.zeros((n, n)); rd = np.zeros((n, n)); W = 0.0
        for i in range(walkers.nw):
            if walkers.w[i] <= 0:
                continue
            try:
                Gu = _green(self.up, walkers.phi_up[i]).T   # <c^+_a c_b>
                Gd = _green(self.dn, walkers.phi_dn[i]).T
            except np.linalg.LinAlgError:
                continue
            w = walkers.w[i]
            ru += w * Gu; rd += w * Gd; W += w
        if W <= 0:
            return
        ru /= W; rd /= W
        ru = 0.5 * (ru + ru.T); rd = 0.5 * (rd + rd.T)   # mixed RDM not Hermitian
        # damp toward the current trial's occupied-orbital projector
        Pu = self.up @ self.up.T; Pd = self.dn @ self.dn.T
        ru = self.mix * Pu + (1.0 - self.mix) * ru
        rd = self.mix * Pd + (1.0 - self.mix) * rd
        # natural orbitals = highest-occupation eigenvectors of the 1-RDM
        _, Vu = np.linalg.eigh(ru); _, Vd = np.linalg.eigh(rd)
        self.up = Vu[:, -self.nup:].copy()
        self.dn = Vd[:, -self.ndn:].copy()
        self._refresh_olp(walkers)            # importance function changed

    def _update_sample1(self, walkers):
        """Stochastic single-walker trial: draw one walker j with probability
        proportional to its weight and use its (orthonormalized) determinant as the
        new trial -- i.e. the constraint/sign boundary is redefined to a randomly
        chosen member of the current generation. A single determinant, so all the
        single-trial machinery applies unchanged; high variance (cf. the smooth
        natural-orbital average)."""
        w = np.clip(walkers.w, 0.0, None)
        if w.sum() <= 0:
            return
        j = int(self.rng.choice(walkers.nw, p=w / w.sum()))
        qu, _ = np.linalg.qr(walkers.phi_up[j]); self.up = qu.copy()
        qd, _ = np.linalg.qr(walkers.phi_dn[j]); self.dn = qd.copy()
        self._refresh_olp(walkers)


class WalkerEnsemble:
    """The population of walker Slater determinants (phi_up, phi_dn), their real
    weights w, and cached trial overlaps olp. Owns stabilization (QR reortho) and
    comb population control."""

    def __init__(self, trial, nw):
        self.nw = nw
        self.phi_up = np.stack([trial.up.copy() for _ in range(nw)])
        self.phi_dn = np.stack([trial.dn.copy() for _ in range(nw)])
        self.w = np.ones(nw)
        self.olp = np.array([trial.overlap(self.phi_up[i], self.phi_dn[i])
                             for i in range(nw)])

    def reorthogonalize(self, trial):
        for i in range(self.nw):
            if self.w[i] <= 0:
                continue
            qu, _ = np.linalg.qr(self.phi_up[i]); self.phi_up[i] = qu
            qd, _ = np.linalg.qr(self.phi_dn[i]); self.phi_dn[i] = qd
        self.olp = np.array([trial.overlap(self.phi_up[i], self.phi_dn[i])
                             for i in range(self.nw)])

    def pop_control(self, rng):
        w = np.clip(self.w, 0, None)
        if w.sum() <= 0:
            return
        w *= self.nw / w.sum()
        # comb resampling
        idx = []
        c = rng.random() / self.nw
        cum = np.cumsum(w) / w.sum()
        for j in range(self.nw):
            target = c + j / self.nw
            idx.append(int(np.searchsorted(cum, target)))
        idx = np.clip(idx, 0, self.nw - 1)
        self.phi_up = self.phi_up[idx].copy()
        self.phi_dn = self.phi_dn[idx].copy()
        self.olp = self.olp[idx].copy()
        self.w = np.ones(self.nw)


class Propagator:
    """One CPMC propagation step: half-K, every density-density term with an
    importance-sampled discrete HS field (constrained-path clip), half-K.
    `step_record` additionally stores the chosen fields for back-propagation."""

    def __init__(self, model, trial, rng):
        self.m = model
        self.trial = trial
        self.rng = rng

    def _term_ratio(self, pu, pd, ov_up, ov_dn, a, sa, b, sb, fa, fb):
        """Overlap ratio after scaling row a of spin-sa det by fa and row b of
        spin-sb det by fb (brute-force determinant ratio)."""
        tr = self.trial
        if sa == sb:
            phi = pu if sa == 0 else pd
            left = tr.up if sa == 0 else tr.dn
            ov0 = ov_up if sa == 0 else ov_dn
            ph = phi.copy(); ph[a, :] *= fa; ph[b, :] *= fb
            return np.linalg.det(left.T @ ph) / ov0
        pha = (pu if sa == 0 else pd).copy(); pha[a, :] *= fa
        lefta = tr.up if sa == 0 else tr.dn
        Ra = np.linalg.det(lefta.T @ pha) / (ov_up if sa == 0 else ov_dn)
        phb = (pu if sb == 0 else pd).copy(); phb[b, :] *= fb
        leftb = tr.up if sb == 0 else tr.dn
        Rb = np.linalg.det(leftb.T @ phb) / (ov_up if sb == 0 else ov_dn)
        return Ra * Rb

    def _step_md_one(self, walkers, i, rec=None):
        """One propagation step for a multi-determinant trial: the importance/
        constraint uses the full (non-factorized) overlap sum_m c_m <up_m|pu><dn_m|pd>,
        recomputed for each candidate field. Slower (loops over the k trial
        determinants) but the single-determinant fast path is left untouched."""
        m, tr = self.m, self.trial
        if walkers.w[i] <= 0:
            if rec is not None: rec[i].append(None)
            return
        pu = m.expK @ walkers.phi_up[i]
        pd = m.expK_dn @ walkers.phi_dn[i]
        ov = tr.md_overlap(pu, pd)
        if ov <= 0:
            walkers.w[i] = 0.0
            if rec is not None: rec[i].append(None)
            return
        choices = []; killed = False
        for (a, sa, b, sb, e1, e2, cf) in m.terms:
            ratios = np.empty(2); cand = []
            for k in (0, 1):
                pu2 = pu.copy(); pd2 = pd.copy()
                (pu2 if sa == 0 else pd2)[a, :] *= e1[k]
                (pu2 if sb == 0 else pd2)[b, :] *= e2[k]
                ov2 = tr.md_overlap(pu2, pd2)
                ratios[k] = cf[k] * max(ov2 / ov, 0.0)
                cand.append((pu2, pd2, ov2))
            tot = ratios.sum()
            if tot <= 0:
                walkers.w[i] = 0.0; killed = True; break
            k = self.rng.choice(2, p=ratios / tot)
            pu, pd, ov = cand[k]
            walkers.w[i] *= 0.5 * tot
            choices.append(k)
        if killed:
            if rec is not None: rec[i].append(None)
            return
        pu = m.expK @ pu; pd = m.expK_dn @ pd
        walkers.phi_up[i] = pu; walkers.phi_dn[i] = pd
        walkers.olp[i] = tr.md_overlap(pu, pd)
        if rec is not None: rec[i].append(choices)

    def step(self, walkers):
        m, tr = self.m, self.trial
        if tr.multi:
            for i in range(walkers.nw):
                self._step_md_one(walkers, i)
            return
        for i in range(walkers.nw):
            if walkers.w[i] <= 0:
                continue
            pu = m.expK @ walkers.phi_up[i]
            pd = m.expK_dn @ walkers.phi_dn[i]
            ov_up = _ov_spin(tr.up, pu)
            ov_dn = _ov_spin(tr.dn, pd)
            if ov_up == 0 or ov_dn == 0:
                walkers.w[i] = 0.0; continue
            killed = False
            for (a, sa, b, sb, e1, e2, cf) in m.terms:
                ratios = np.empty(2)
                for k in (0, 1):                  # HS field s=+1 (k=0), -1 (k=1)
                    R = self._term_ratio(pu, pd, ov_up, ov_dn, a, sa, b, sb, e1[k], e2[k])
                    ratios[k] = cf[k] * max(R, 0.0)   # constrained path: clip R<=0
                tot = ratios.sum()
                if tot <= 0:
                    walkers.w[i] = 0.0; killed = True; break
                k = self.rng.choice(2, p=ratios / tot)
                if sa == 0: pu[a, :] *= e1[k]
                else:       pd[a, :] *= e1[k]
                if sb == 0: pu[b, :] *= e2[k]
                else:       pd[b, :] *= e2[k]
                ov_up = _ov_spin(tr.up, pu)
                ov_dn = _ov_spin(tr.dn, pd)
                walkers.w[i] *= 0.5 * tot
            if killed:
                continue
            pu = m.expK @ pu
            pd = m.expK_dn @ pd
            walkers.phi_up[i] = pu; walkers.phi_dn[i] = pd
            walkers.olp[i] = tr.overlap(pu, pd)

    def step_record(self, walkers, rec):
        """One propagation step, recording per-walker the HS field chosen for each
        term (so the path can be replayed for back-propagation). rec[i] gets a list
        of field indices (one per term, in model.terms order), or None if the
        walker died this step."""
        m, tr = self.m, self.trial
        if tr.multi:
            for i in range(walkers.nw):
                self._step_md_one(walkers, i, rec)
            return
        for i in range(walkers.nw):
            if walkers.w[i] <= 0:
                rec[i].append(None); continue
            pu = m.expK @ walkers.phi_up[i]
            pd = m.expK_dn @ walkers.phi_dn[i]
            ov_up = _ov_spin(tr.up, pu)
            ov_dn = _ov_spin(tr.dn, pd)
            if ov_up == 0 or ov_dn == 0:
                walkers.w[i] = 0.0; rec[i].append(None); continue
            choices = []; killed = False
            for (a, sa, b, sb, e1, e2, cf) in m.terms:
                ratios = np.empty(2)
                for k in (0, 1):
                    R = self._term_ratio(pu, pd, ov_up, ov_dn, a, sa, b, sb, e1[k], e2[k])
                    ratios[k] = cf[k] * max(R, 0.0)
                tot = ratios.sum()
                if tot <= 0:
                    walkers.w[i] = 0.0; killed = True; break
                k = self.rng.choice(2, p=ratios / tot)
                if sa == 0: pu[a, :] *= e1[k]
                else:       pd[a, :] *= e1[k]
                if sb == 0: pu[b, :] *= e2[k]
                else:       pd[b, :] *= e2[k]
                ov_up = _ov_spin(tr.up, pu)
                ov_dn = _ov_spin(tr.dn, pd)
                walkers.w[i] *= 0.5 * tot
                choices.append(k)
            if killed:
                rec[i].append(None); continue
            pu = m.expK @ pu
            pd = m.expK_dn @ pd
            walkers.phi_up[i] = pu; walkers.phi_dn[i] = pd
            rec[i].append(choices)


class Estimators:
    """Observable estimators given walker Green's functions: mixed-estimator
    energy, and back-propagated energy + equal-time correlations."""

    def __init__(self, model, trial):
        self.m = model
        self.trial = trial

    def _two_body(self, G):
        """sum over pot_terms of V <n_a n_b> via Wick (G = (Gu, Gd))."""
        ev = 0.0
        for (a, sa, b, sb, V) in self.m.pot_terms:
            if sa != sb:                       # different spin: no exchange
                ev += V * G[sa][a, a] * G[sb][b, b]
            else:                              # same spin: Wick exchange
                Gs = G[sa]
                ev += V * (Gs[a, a] * Gs[b, b] - Gs[a, b] * Gs[b, a])
        return ev

    def _local_energy(self, phi_up, phi_dn):
        """Local energy E_L = <psi_T|H|phi>/<psi_T|phi>. For a multi-determinant
        trial the TWO-BODY term is quadratic in G, so it must be the t_m-weighted
        average of the PER-DETERMINANT Wick energies (NOT the Wick form of the
        averaged Green's function): E_L = sum_m t_m E^(m) / sum_m t_m,
        t_m = c_m <up_m|phi><dn_m|phi>."""
        m, tr = self.m, self.trial
        if not tr.multi:
            Gu = _green(tr.up, phi_up); Gd = _green(tr.dn, phi_dn)
            return np.sum(m.K * Gu.T) + np.sum(m.K_dn * Gd.T) + self._two_body((Gu, Gd))
        numw = 0.0; denw = 0.0
        for mm in range(len(tr.coef)):
            ou = _ov_spin(tr.ups[mm], phi_up); od = _ov_spin(tr.dns[mm], phi_dn)
            t = tr.coef[mm] * ou * od
            if t == 0.0:
                continue
            Gu = _green(tr.ups[mm], phi_up); Gd = _green(tr.dns[mm], phi_dn)
            em = np.sum(m.K * Gu.T) + np.sum(m.K_dn * Gd.T) + self._two_body((Gu, Gd))
            numw += t * em; denw += t
        return numw / denw if denw != 0 else float("nan")

    def energy(self, walkers):
        num = 0.0; den = 0.0
        for i in range(walkers.nw):
            if walkers.w[i] <= 0:
                continue
            el = self._local_energy(walkers.phi_up[i], walkers.phi_dn[i])
            if not np.isfinite(el):
                continue
            num += walkers.w[i] * el
            den += walkers.w[i]
        return num / den if den > 0 else float("nan")

    def bp_bra(self, rec_i, bp):
        """Back-propagated bra: apply the recorded B_step operators to psiT in
        reverse time order (B^dagger = B for our Hermitian symmetric/diagonal ops)."""
        m, tr = self.m, self.trial
        Lu = tr.up.copy(); Ld = tr.dn.copy()
        for s in range(bp - 1, -1, -1):
            ch = rec_i[s]
            if ch is None:
                return None, None                      # walker died in the window
            Lu = m.expK @ Lu; Ld = m.expK_dn @ Ld
            for (a, sa, b, sb, e1, e2, cf), k in zip(m.terms, ch):
                if sa == 0: Lu[a, :] *= e1[k]
                else:       Ld[a, :] *= e1[k]
                if sb == 0: Lu[b, :] *= e2[k]
                else:       Ld[b, :] *= e2[k]
            Lu = m.expK @ Lu; Ld = m.expK_dn @ Ld
        return Lu, Ld

    def energy_bp(self, walkers, ket_up, ket_dn, rec, bp):
        m = self.m
        num = den = 0.0
        for i in range(walkers.nw):
            if walkers.w[i] <= 0:
                continue
            Lu, Ld = self.bp_bra(rec[i], bp)
            if Lu is None:
                continue
            try:
                Gu = _green(Lu, ket_up[i]); Gd = _green(Ld, ket_dn[i])
            except np.linalg.LinAlgError:
                continue
            ek = np.sum(m.K * Gu.T) + np.sum(m.K_dn * Gd.T)
            ev = self._two_body((Gu, Gd))
            num += walkers.w[i] * (ek + ev); den += walkers.w[i]
        return num / den if den > 0 else float("nan")

    def corr_block(self, walkers, ket_up, ket_dn, rec, bp):
        """Per-block weighted sums of the back-propagated equal-time observables:
        G^s_ij = <c^+_{i,s} c_{j,s}>, charge <n_i n_j>, spin <S^z_i S^z_j>
        (Wick from each walker's BP Green's function). Returns (Gu, Gd, nn, sz, W)."""
        n = self.m.n
        Gu_s = np.zeros((n, n)); Gd_s = np.zeros((n, n))
        nn_s = np.zeros((n, n)); sz_s = np.zeros((n, n)); W = 0.0
        for i in range(walkers.nw):
            if walkers.w[i] <= 0:
                continue
            Lu, Ld = self.bp_bra(rec[i], bp)
            if Lu is None:
                continue
            try:
                Gu = _green(Lu, ket_up[i]).T   # Gu[i,j] = <c^+_i c_j>
                Gd = _green(Ld, ket_dn[i]).T
            except np.linalg.LinAlgError:
                continue
            w = walkers.w[i]
            nu = np.diag(Gu); nd = np.diag(Gd)
            exu = np.diag(nu) - Gu * Gu.T            # delta_ij n_i - G_ij G_ji
            exd = np.diag(nd) - Gd * Gd.T
            nn = np.outer(nu + nd, nu + nd) + exu + exd
            sz = 0.25 * (np.outer(nu - nd, nu - nd) + exu + exd)
            Gu_s += w * Gu; Gd_s += w * Gd; nn_s += w * nn; sz_s += w * sz; W += w
        return Gu_s, Gd_s, nn_s, sz_s, W

    def pairmag_block(self, walkers, ket_up, ket_dn, rec, bp, Bs, Bd, phase):
        """Per-block weighted back-propagated pairing + magnetic observables:
          Ps/Pd : singlet pairing structure factor S_alpha = sum_{m,n} P_alpha(m,n)
                  with P_alpha(m,n) = Gu[m,n] (B_alpha Gd B_alpha^T)[m,n]
                  (alpha = s, d_{x^2-y^2}; B_alpha = bond form-factor matrix);
          Sq    : spin structure factor (1/N) phase^T (S^z S^z) phase  (q from phase);
          mom   : weighted sum of the local-moment vector m_i^z = n_iu - n_id
                  (caller squares the average -> sum_i <m_i^z>^2).
        Wick-contracted from each walker's BP Green's function. Also accumulates the
        weighted Green's functions Gu_s, Gd_s so the caller can form the DISCONNECTED
        pairing (from the averaged G) and hence the CONNECTED/VERTEX pairing
        P_vertex = P_full - P_disc. Returns (Ps, Pd, Sq, mom_vec, Gu_s, Gd_s, W)."""
        n = self.m.n
        lxy = self.m.lx * self.m.ly
        gAM = np.where(np.arange(n) < lxy, 1.0, -1.0)   # orbital-staggered sign (AM order)
        Ps = Pd = Sq = Sam = W = 0.0
        mom = np.zeros(n); Gu_s = np.zeros((n, n)); Gd_s = np.zeros((n, n))
        for i in range(walkers.nw):
            if walkers.w[i] <= 0:
                continue
            Lu, Ld = self.bp_bra(rec[i], bp)
            if Lu is None:
                continue
            try:
                Gu = _green(Lu, ket_up[i]).T   # Gu[i,j] = <c^+_i c_j>
                Gd = _green(Ld, ket_dn[i]).T
            except np.linalg.LinAlgError:
                continue
            w = walkers.w[i]
            # pairing (FULL): P_alpha(m,n) = Gu[m,n] * (B_alpha Gd B_alpha^T)[m,n]
            Ps += w * float((Gu * (Bs @ Gd @ Bs.T)).sum())
            Pd += w * float((Gu * (Bd @ Gd @ Bd.T)).sum())
            Gu_s += w * Gu; Gd_s += w * Gd
            # spin S^z S^z (Wick), then structure factor at the phase's q
            nu = np.diag(Gu); nd = np.diag(Gd)
            exu = np.diag(nu) - Gu * Gu.T
            exd = np.diag(nd) - Gd * Gd.T
            sz = 0.25 * (np.outer(nu - nd, nu - nd) + exu + exd)
            Sq += w * float(phase @ sz @ phase) / n
            Sam += w * float(gAM @ sz @ gAM) / n        # altermagnetic (orbital-staggered) SF
            mom += w * (nu - nd)
            W += w
        return Ps, Pd, Sq, Sam, mom, Gu_s, Gd_s, W

    def _step_bmats(self, ch):
        """Reconstruct the per-spin one-body propagators b^sigma = e^{-dtK/2} D^sigma
        e^{-dtK/2} (and their inverses) for one recorded step, from the chosen HS
        field index list `ch` (in model.terms order). D^sigma is the diagonal of the
        accumulated HS row factors for that spin."""
        m = self.m
        n = m.n
        d_up = np.ones(n); d_dn = np.ones(n)
        for (a, sa, b, sb, e1, e2, cf), k in zip(m.terms, ch):
            (d_up if sa == 0 else d_dn)[a] *= e1[k]
            (d_up if sb == 0 else d_dn)[b] *= e2[k]
        bu = m.expK @ (d_up[:, None] * m.expK)
        bd = m.expK_dn @ (d_dn[:, None] * m.expK_dn)
        bui = m.expK_inv @ ((1.0 / d_up)[:, None] * m.expK_inv)
        bdi = m.expK_dn_inv @ ((1.0 / d_dn)[:, None] * m.expK_dn_inv)
        return bu, bd, bui, bdi

    def chi_block(self, walkers, ket_up, ket_dn, rec, bp, phase):
        """Imaginary-time-displaced staggered spin correlation C(tau_l) =
        <O(tau_l) O(0)>, O = sum_i phase_i S^z_i, for l = 0..bp. The ket sits at
        slice 0, the recorded fields give the forward propagators B_(l), and the
        back-propagated bra fixes the ground-state equal-time GF g^s_ij=<c^+_i c_j>.
        Per spin: P^s(tau)=B_(l)(I-g^sT) [particle], H^s(tau)=g^s B_(l)^{-1} [hole],
        g^s(tau)=B_(l) g^s B_(l)^{-1} [equal-time at tau]. Returns (Csum[0..bp], W)."""
        L = bp + 1
        Csum = np.zeros(L); W = 0.0
        for i in range(walkers.nw):
            if walkers.w[i] <= 0:
                continue
            chs = rec[i]
            if any(c is None for c in chs):
                continue
            Lu, Ld = self.bp_bra(chs, bp)
            if Lu is None:
                continue
            try:
                gu = _green(Lu, ket_up[i]).T   # g_ij = <c^+_i c_j>
                gd = _green(Ld, ket_dn[i]).T
            except np.linalg.LinAlgError:
                continue
            w = walkers.w[i]
            n = self.m.n
            Bu = np.eye(n); Bd = np.eye(n); Bui = np.eye(n); Bdi = np.eye(n)
            ImguT = np.eye(n) - gu.T; ImgdT = np.eye(n) - gd.T
            O0 = 0.5 * (phase @ np.diag(gu) - phase @ np.diag(gd))   # O(0)
            for l in range(L):
                if l > 0:
                    bu, bd, bui, bdi = self._step_bmats(chs[l - 1])
                    Bu = bu @ Bu; Bd = bd @ Bd
                    Bui = Bui @ bui; Bdi = Bdi @ bdi
                # equal-time GF at slice l:  g^s(tau) = B^{-T} g B^{T}
                # (c_i(tau)=B c, c_i^dag(tau)=B^{-T} c^dag for B=e^{-tau h})
                nu = np.diag(Bui.T @ gu @ Bu.T); nd = np.diag(Bdi.T @ gd @ Bd.T)
                Otau = 0.5 * (phase @ nu - phase @ nd)
                disc = Otau * O0
                # particle P^s=B(I-g^T)=<c_i(t)c_j^dag>, hole H^s=B^{-T}g=<c_i^dag(t)c_j>
                # same-spin connected: (1/4) sum_s phase^T (H^s o P^s) phase
                Hu = Bui.T @ gu; Pu = Bu @ ImguT
                Hd = Bdi.T @ gd; Pd = Bd @ ImgdT
                conn = 0.25 * (phase @ (Hu * Pu) @ phase + phase @ (Hd * Pd) @ phase)
                Csum[l] += w * (disc + conn)
            W += w
        return Csum, W

    def chi_spin_block(self, walkers, ket_up, ket_dn, rec, bp):
        """Imaginary-time-displaced spin correlation MATRIX for ALL momenta at once:
        M_l[i,j] = <S^z_i(tau_l) S^z_j(0)> for l = 0..bp. This is the exact matrix
        generalization of chi_block (which contracts with one hard-coded phase
        vector): per walker and slice, with m_i = n_iu - n_id,
            M_l[i,j] = 1/4 m_i(tau_l) m_j(0)  +  1/4 sum_s (H^s o P^s)[i,j],
        using the same particle/hole time-displaced GFs P^s(tau)=B_(l)(I-g^sT),
        H^s(tau)=B_(l)^{-T} g^s and equal-time g^s(tau)=B^{-T} g B^{T} as chi_block
        (chi_block's scalar == phase^T M_l phase, bit-checked in the validator).
        Also the TRANSVERSE channel (single cross-spin Wick pairing):
            Mpm_l[i,j] = <S+_i(tau_l) S-_j(0)> = H^up_{ij}(tau_l) P^dn_{ij}(tau_l)
        (SU(2): chi_pm = 2 chi_zz). The caller Fourier-reduces M_l -> C_q(tau_l).
        Returns (Msum[L,n,n], Psum[L,n,n], W)."""
        L = bp + 1; n = self.m.n
        Msum = np.zeros((L, n, n)); Psum = np.zeros((L, n, n)); W = 0.0
        for i in range(walkers.nw):
            if walkers.w[i] <= 0:
                continue
            chs = rec[i]
            if any(c is None for c in chs):
                continue
            Lu, Ld = self.bp_bra(chs, bp)
            if Lu is None:
                continue
            try:
                gu = _green(Lu, ket_up[i]).T   # g_ij = <c^+_i c_j>
                gd = _green(Ld, ket_dn[i]).T
            except np.linalg.LinAlgError:
                continue
            w = walkers.w[i]
            Bu = np.eye(n); Bd = np.eye(n); Bui = np.eye(n); Bdi = np.eye(n)
            ImguT = np.eye(n) - gu.T; ImgdT = np.eye(n) - gd.T
            m0 = 0.5 * (np.diag(gu) - np.diag(gd))              # m_j(0)/2
            for l in range(L):
                if l > 0:
                    bu, bd, bui, bdi = self._step_bmats(chs[l - 1])
                    Bu = bu @ Bu; Bd = bd @ Bd
                    Bui = Bui @ bui; Bdi = Bdi @ bdi
                nu = np.diag(Bui.T @ gu @ Bu.T); nd = np.diag(Bdi.T @ gd @ Bd.T)
                ml = 0.5 * (nu - nd)                            # m_i(tau_l)/2
                Hu = Bui.T @ gu; Pu = Bu @ ImguT                # hole/particle GFs
                Hd = Bdi.T @ gd; Pd = Bd @ ImgdT
                Msum[l] += w * (np.outer(ml, m0) + 0.25 * (Hu * Pu + Hd * Pd))
                Psum[l] += w * (Hu * Pd)                        # <c^+_u(t)c_u><c_d(t)c^+_d>
            W += w
        return Msum, Psum, W

    def chid_block(self, walkers, ket_up, ket_dn, rec, bp, Fs, Fd):
        """Imaginary-time-displaced singlet PAIRING correlation
        C_a(tau_l) = <Delta_a(tau_l) Delta_a^dag(0)>, a = s, d_{x^2-y^2}, for the
        pair operator Delta_a^dag = sum_m sum_delta f_a(delta) c^+_{m up} c^+_{m+d dn}.
        By Wick (spin-separated) C_a(tau) = sum_{m,n} P^up(tau)_{mn} (F_a P^dn(tau)
        F_a^T)_{mn}, with the time-displaced PARTICLE GF P^s(tau)=B_(l)(I-g^sT) and
        F_a the bond form-factor matrix. This is the FULL (bubble + vertex)
        susceptibility kernel. To separate the VERTEX the caller also needs the
        BUBBLE built from the ensemble-AVERAGED time-displaced GFs, so we accumulate
        per-tau the weighted P^up(tau), P^dn(tau) (Pu_s, Pd_s). Returns
        (Cs, Cd, Pu_s, Pd_s, W) with Cs/Cd shape (L,), Pu_s/Pd_s shape (L,n,n)."""
        L = bp + 1; n = self.m.n
        Cs = np.zeros(L); Cd = np.zeros(L); W = 0.0
        Pu_s = np.zeros((L, n, n)); Pd_s = np.zeros((L, n, n))
        for i in range(walkers.nw):
            if walkers.w[i] <= 0:
                continue
            chs = rec[i]
            if any(c is None for c in chs):
                continue
            Lu, Ld = self.bp_bra(chs, bp)
            if Lu is None:
                continue
            try:
                gu = _green(Lu, ket_up[i]).T; gd = _green(Ld, ket_dn[i]).T
            except np.linalg.LinAlgError:
                continue
            w = walkers.w[i]
            Bu = np.eye(n); Bd = np.eye(n)
            ImguT = np.eye(n) - gu.T; ImgdT = np.eye(n) - gd.T
            for l in range(L):
                if l > 0:
                    bu, bd, _, _ = self._step_bmats(chs[l - 1])
                    Bu = bu @ Bu; Bd = bd @ Bd
                Pu = Bu @ ImguT; Pd = Bd @ ImgdT     # time-displaced particle GFs
                Cs[l] += w * float(np.sum(Pu * (Fs @ Pd @ Fs.T)))
                Cd[l] += w * float(np.sum(Pu * (Fd @ Pd @ Fd.T)))
                Pu_s[l] += w * Pu; Pd_s[l] += w * Pd
            W += w
        return Cs, Cd, Pu_s, Pd_s, W


class CPMC:
    """Orchestrator: builds the model, trial, walker ensemble, propagator and
    estimators, and drives equilibration + measurement. Public surface (run,
    run_bp, run_bp_obs) and constructor are unchanged across the Phase 6 refactor."""

    def __init__(self, lx=4, ly=4, nup=1, ndn=1, t=1.0, U=3.0, dt=0.01,
                 nwalkers=400, seed=1, K=None, uxy=0.0, v=0.0,
                 trial="fixed", trial_every=10, trial_mix=0.5, trial_k=8, K_dn=None):
        # K: optional custom hopping matrix (n x n). If None, single-band square
        # lattice. For the two-orbital altermagnet pass build_hopping(...); the
        # density-density machinery is dimension-agnostic. K_dn: optional separate
        # down-spin hopping (SPIN-DEPENDENT model, e.g. the altermagnet tam model);
        # None => spin-independent (K_dn=K).
        # trial: "fixed" (free-electron GS) or "adaptive" (self-consistent, WS3);
        #        trial_every = update cadence in steps; trial_mix = damping.
        self.lx, self.ly = lx, ly
        Kmat = square_hopping(lx, ly, t) if K is None else np.asarray(K, float)
        self.nup, self.ndn, self.U, self.dt = nup, ndn, U, dt
        self.uxy, self.v = uxy, v
        self.nw = nwalkers
        self.trial_every = trial_every
        self.rng = np.random.default_rng(seed)

        self.model = LatticeModel(Kmat, lx, ly, dt, U, uxy, v, K_dn=K_dn)
        self.n = self.model.n
        # trial = free-electron ground state (lowest orbitals of each spin's K)
        vu = self.model.eigvecs; vd = self.model.eigvecs_dn
        self.trial = TrialWF(vu[:, :nup].copy(), vd[:, :ndn].copy(),
                             mode=trial, mix=trial_mix, rng=self.rng, k=trial_k)
        self.walkers = WalkerEnsemble(self.trial, self.nw)
        self.prop = Propagator(self.model, self.trial, self.rng)
        self.est = Estimators(self.model, self.trial)

    def _maybe_update_trial(self, it):
        if self.trial.mode != "fixed" and (it + 1) % self.trial_every == 0:
            self.trial.update(self.walkers)

    # ---- convenience accessors (state lives in the sub-objects) ----
    @property
    def K(self): return self.model.K
    @property
    def phi_up(self): return self.walkers.phi_up
    @property
    def phi_dn(self): return self.walkers.phi_dn
    @property
    def w(self): return self.walkers.w

    def step(self):
        self.prop.step(self.walkers)

    def reorthogonalize(self):
        self.walkers.reorthogonalize(self.trial)

    def pop_control(self):
        self.walkers.pop_control(self.rng)

    def energy(self):
        return self.est.energy(self.walkers)

    def run(self, nequil=200, nmeas=400, ortho=10, pc=10, meas_every=5):
        for it in range(nequil):
            self.step()
            if (it + 1) % ortho == 0:
                self.reorthogonalize()
            if (it + 1) % pc == 0:
                self.pop_control()
            self._maybe_update_trial(it)
        es = []
        for it in range(nmeas):
            self.step()
            if (it + 1) % ortho == 0:
                self.reorthogonalize()
            if (it + 1) % pc == 0:
                self.pop_control()
            self._maybe_update_trial(it)
            if (it + 1) % meas_every == 0:
                es.append(self.energy())
        es = np.array(es)
        return es.mean(), es.std() / np.sqrt(len(es))

    def run_bp(self, nequil=150, nblocks=40, bp=15, ortho=10, pc=10):
        """Ground-state CPMC with back-propagated energy (less biased than the
        mixed estimator for operators not commuting with H, e.g. uxy/v)."""
        for it in range(nequil):
            self.step()
            if (it + 1) % ortho == 0: self.reorthogonalize()
            if (it + 1) % pc == 0: self.pop_control()
            self._maybe_update_trial(it)
        es = []
        for blk in range(nblocks):
            self.reorthogonalize()
            ket_up = self.walkers.phi_up.copy(); ket_dn = self.walkers.phi_dn.copy()
            rec = [[] for _ in range(self.nw)]
            for _ in range(bp):
                self.prop.step_record(self.walkers, rec)
            es.append(self.est.energy_bp(self.walkers, ket_up, ket_dn, rec, bp))
            self.pop_control()
            # adaptive/sampled trial updated between blocks (frozen across BP window)
            if self.trial.mode != "fixed":
                self.trial.update(self.walkers)
        es = np.array([e for e in es if np.isfinite(e)])
        return es.mean(), es.std() / np.sqrt(len(es))

    def run_bp_obs(self, nequil=150, nblocks=30, bp=12, ortho=10, pc=10):
        """Back-propagated equal-time observables: G^s_ij = <c^+_{i,s} c_{j,s}>,
        charge <n_i n_j>, spin <S^z_i S^z_j> (Wick from each walker's BP GF)."""
        n = self.n
        for it in range(nequil):
            self.step()
            if (it + 1) % ortho == 0: self.reorthogonalize()
            if (it + 1) % pc == 0: self.pop_control()
        Gu_s = np.zeros((n, n)); Gd_s = np.zeros((n, n))
        nn_s = np.zeros((n, n)); sz_s = np.zeros((n, n)); W = 0.0
        for blk in range(nblocks):
            self.reorthogonalize()
            ket_up = self.walkers.phi_up.copy(); ket_dn = self.walkers.phi_dn.copy()
            rec = [[] for _ in range(self.nw)]
            for _ in range(bp):
                self.prop.step_record(self.walkers, rec)
            bGu, bGd, bnn, bsz, bW = self.est.corr_block(self.walkers, ket_up, ket_dn, rec, bp)
            Gu_s += bGu; Gd_s += bGd; nn_s += bnn; sz_s += bsz; W += bW
            self.pop_control()
        return {"nsites": n, "green_up": (Gu_s / W).tolist(), "green_dn": (Gd_s / W).tolist(),
                "nn": (nn_s / W).tolist(), "szsz": (sz_s / W).tolist()}

    def _square_bond_factors(self):
        """For the single-band square lattice (index i = x*ly + y): the bond
        form-factor matrices B_s (s-wave, all +1) and B_d (d_{x^2-y^2}: +1 on x
        bonds, -1 on y bonds) with B[m, neighbour] = f, and the (pi,pi) phase."""
        lx, ly, n = self.lx, self.ly, self.n
        assert n == lx * ly, "run_bp_pairmag is for the single-band square lattice"
        Bs = np.zeros((n, n)); Bd = np.zeros((n, n)); phase = np.zeros(n)
        for x in range(lx):
            for y in range(ly):
                m = x * ly + y
                phase[m] = (-1.0) ** (x + y)
                nbrs = [(((x + 1) % lx) * ly + y, +1.0),   # +x  (d: +1)
                        (((x - 1) % lx) * ly + y, +1.0),   # -x  (d: +1)
                        (x * ly + (y + 1) % ly, -1.0),     # +y  (d: -1)
                        (x * ly + (y - 1) % ly, -1.0)]     # -y  (d: -1)
                for (j, fd) in nbrs:
                    Bs[m, j] += 1.0
                    Bd[m, j] += fd
        return Bs, Bd, phase

    def _bond_factors(self):
        """Bond form-factor matrices for the single-band square lattice OR the
        two-orbital altermagnet (n = norb * lx*ly, index = orb*lxy + x*ly + y).
        Intra-orbital nearest-neighbour s-wave (all +1) and d_{x^2-y^2} (+1 on x
        bonds, -1 on y bonds), plus the (pi,pi) staggered phase. For norb=1 this is
        identical to _square_bond_factors."""
        lx, ly = self.lx, self.ly
        lxy = lx * ly
        n = self.n
        norb = n // lxy
        assert norb * lxy == n, "n must be norb * lx * ly"
        Fs = np.zeros((n, n)); Fd = np.zeros((n, n)); phase = np.zeros(n)
        def idx(x, y, orb): return orb * lxy + (x % lx) * ly + (y % ly)
        for orb in range(norb):
            for x in range(lx):
                for y in range(ly):
                    m = idx(x, y, orb)
                    phase[m] = (-1.0) ** (x + y)
                    nbrs = [(idx(x + 1, y, orb), +1.0), (idx(x - 1, y, orb), +1.0),
                            (idx(x, y + 1, orb), -1.0), (idx(x, y - 1, orb), -1.0)]
                    for (j, fd) in nbrs:
                        Fs[m, j] += 1.0
                        Fd[m, j] += fd
        return Fs, Fd, phase

    def run_bp_pairmag(self, nequil=150, nblocks=30, bp=12, ortho=10, pc=10):
        """Back-propagated singlet pairing (s + d_{x^2-y^2}) and magnetic observables
        for the single-band square lattice. Reports the FULL pairing structure factor
        P_full = <Delta Delta^dag>, the DISCONNECTED part (from the averaged Green's
        functions), and the CONNECTED/VERTEX part P_vertex = P_full - P_disc (the
        quantity used in the manuscript), plus S(pi,pi) and local moment^2. Per-block
        error bars; returns {value, error} dicts."""
        Bs, Bd, phase = self._bond_factors()
        for it in range(nequil):
            self.step()
            if (it + 1) % ortho == 0: self.reorthogonalize()
            if (it + 1) % pc == 0: self.pop_control()
            self._maybe_update_trial(it)         # let adaptive/CASCI improve the trial
        keys = ["pair_swave", "pair_dwave", "pair_swave_vertex", "pair_dwave_vertex",
                "Sq_pipi", "S_AM", "moment2"]
        bvals = {k: [] for k in keys}
        for blk in range(nblocks):
            self.reorthogonalize()
            ket_up = self.walkers.phi_up.copy(); ket_dn = self.walkers.phi_dn.copy()
            rec = [[] for _ in range(self.nw)]
            for _ in range(bp):
                self.prop.step_record(self.walkers, rec)
            bPs, bPd, bSq, bSam, bmom, bGu, bGd, bW = self.est.pairmag_block(
                self.walkers, ket_up, ket_dn, rec, bp, Bs, Bd, phase)
            if bW > 0:
                Gu = bGu / bW; Gd = bGd / bW                     # averaged G -> disc
                disc_s = float((Gu * (Bs @ Gd @ Bs.T)).sum())
                disc_d = float((Gu * (Bd @ Gd @ Bd.T)).sum())
                bvals["pair_swave"].append(bPs / bW)
                bvals["pair_dwave"].append(bPd / bW)
                bvals["pair_swave_vertex"].append(bPs / bW - disc_s)
                bvals["pair_dwave_vertex"].append(bPd / bW - disc_d)
                bvals["Sq_pipi"].append(bSq / bW)
                bvals["S_AM"].append(bSam / bW)
                bvals["moment2"].append(float(np.sum((bmom / bW) ** 2)))
            self.pop_control()
        out = {"nsites": self.n}
        for k, vs in bvals.items():
            a = np.array(vs)
            out[k] = {"value": float(a.mean()),
                      "error": float(a.std() / np.sqrt(len(a))) if len(a) > 1 else 0.0}
        return out

    def run_bp_chi(self, nequil=150, nblocks=30, bp=20, ortho=10, pc=10):
        """Unequal-time (imaginary-time-displaced) STAGGERED spin susceptibility for
        the single-band square lattice: the correlation C(tau_l)=<O(tau_l) O(0)>,
        O = sum_i (-1)^{x+y} S^z_i, for tau_l = l*dt (l=0..bp), and the static
        susceptibility chi_s = integral_0^inf C(tau) dtau ~ dt * sum_l C(tau_l)
        (trapezoidal). Per-block estimates give error bars. C(tau_0) equals the
        equal-time staggered structure factor N*S(pi,pi) (consistency check)."""
        _, _, phase = self._square_bond_factors()
        for it in range(nequil):
            self.step()
            if (it + 1) % ortho == 0: self.reorthogonalize()
            if (it + 1) % pc == 0: self.pop_control()
        blocks = []
        for blk in range(nblocks):
            self.reorthogonalize()
            ket_up = self.walkers.phi_up.copy(); ket_dn = self.walkers.phi_dn.copy()
            rec = [[] for _ in range(self.nw)]
            for _ in range(bp):
                self.prop.step_record(self.walkers, rec)
            Csum, W = self.est.chi_block(self.walkers, ket_up, ket_dn, rec, bp, phase)
            if W > 0:
                blocks.append(Csum / W)
            self.pop_control()
        B = np.array(blocks)                       # (nblocks, bp+1)
        Ctau = B.mean(axis=0)
        Cerr = B.std(axis=0) / np.sqrt(len(B)) if len(B) > 1 else np.zeros_like(Ctau)
        taus = self.dt * np.arange(bp + 1)
        # chi = trapezoidal integral over tau, per block then averaged (for error)
        chi_blocks = self.dt * (B[:, 1:-1].sum(axis=1) + 0.5 * (B[:, 0] + B[:, -1]))
        return {"nsites": self.n, "bp": bp, "dt": self.dt,
                "taus": taus.tolist(), "Ctau": Ctau.tolist(), "Cerr": Cerr.tolist(),
                "chi_stag": float(chi_blocks.mean()),
                "chi_stag_err": float(chi_blocks.std() / np.sqrt(len(chi_blocks))
                                      if len(chi_blocks) > 1 else 0.0)}

    def run_bp_chi_spin(self, nequil=150, nblocks=30, bp=20, ortho=10, pc=10):
        """Momentum-resolved unequal-time spin correlation for the single-band
        square lattice: C_q(tau_l) = (1/N) sum_ij e^{-iq(ri-rj)} <S^z_i(tau_l) S^z_j(0)>
        on the full (lx,ly) q-grid (reduce_mat FFT convention, identical to the
        DQMC/CP-DQMC grids), and the WINDOWED static susceptibility
        chi_zz(q) = trapezoid_{0..bp*dt} C_q(tau) dtau (per-site; one-sided T=0 Kubo,
        same windowed convention as run_bp_chi). Per-block error bars.
        Consistency: N * C_q(pi,pi) == run_bp_chi's staggered C(tau)."""
        assert self.n == self.lx * self.ly, "run_bp_chi_spin is single-band only"
        from unified_scan import _shift_index
        for it in range(nequil):
            self.step()
            if (it + 1) % ortho == 0: self.reorthogonalize()
            if (it + 1) % pc == 0: self.pop_control()
        blocks = []; blocks_pm = []
        for blk in range(nblocks):
            self.reorthogonalize()
            ket_up = self.walkers.phi_up.copy(); ket_dn = self.walkers.phi_dn.copy()
            rec = [[] for _ in range(self.nw)]
            for _ in range(bp):
                self.prop.step_record(self.walkers, rec)
            Msum, Psum, W = self.est.chi_spin_block(self.walkers, ket_up, ket_dn, rec, bp)
            if W > 0:
                blocks.append(Msum / W); blocks_pm.append(Psum / W)
            self.pop_control()
        shift = _shift_index(self.lx, self.ly)
        rows = np.arange(self.n)
        def _pq(M):     # M[i,j] -> S(R)=sum_m M[m,m+R] -> P(q)=FFT2(S)  (== reduce_mat's Pq)
            Sg = np.array([[M[rows, shift[dx, dy]].sum() for dy in range(self.ly)]
                           for dx in range(self.lx)])
            return np.real(np.fft.fft2(Sg))
        # per block, per tau: full q-grid C_q(tau_l) = FFT2(S(R))/N
        Cq = np.array([[_pq(Mb[l]) / self.n
                        for l in range(bp + 1)] for Mb in blocks])   # (nb, L, lx, ly)
        nb = len(blocks)
        Ctau_q = Cq.mean(axis=0)
        Cerr_q = Cq.std(axis=0) / np.sqrt(nb) if nb > 1 else np.zeros_like(Ctau_q)
        # windowed chi(q): trapezoid per block, then block stats
        chi_b = self.dt * (Cq[:, 1:-1].sum(axis=1) + 0.5 * (Cq[:, 0] + Cq[:, -1]))
        chi_q = chi_b.mean(axis=0)
        chi_q_err = chi_b.std(axis=0) / np.sqrt(nb) if nb > 1 else np.zeros_like(chi_q)
        # transverse channel: same reduction/window on the S+S- matrices
        Pq = np.array([[_pq(Mb[l]) / self.n
                        for l in range(bp + 1)] for Mb in blocks_pm])
        Ptau_q = Pq.mean(axis=0)
        Perr_q = Pq.std(axis=0) / np.sqrt(nb) if nb > 1 else np.zeros_like(Ptau_q)
        pm_b = self.dt * (Pq[:, 1:-1].sum(axis=1) + 0.5 * (Pq[:, 0] + Pq[:, -1]))
        chi_pm_q = pm_b.mean(axis=0)
        chi_pm_err = pm_b.std(axis=0) / np.sqrt(nb) if nb > 1 else np.zeros_like(chi_pm_q)
        kpi = (self.lx // 2, self.ly // 2)                     # (pi,pi) grid index
        return {"nsites": self.n, "bp": bp, "dt": self.dt,
                "taus": (self.dt * np.arange(bp + 1)).tolist(),
                "Ctau_q": Ctau_q.tolist(), "Cerr_q": Cerr_q.tolist(),
                "chi_q": chi_q.tolist(), "chi_q_err": chi_q_err.tolist(),
                "Ctau_pm_q": Ptau_q.tolist(), "Cerr_pm_q": Perr_q.tolist(),
                "chi_pm_q": chi_pm_q.tolist(), "chi_pm_q_err": chi_pm_err.tolist(),
                "chi_q0": float(chi_q[0, 0]),
                "chi_pipi": float(chi_q[kpi]), "chi_pipi_err": float(chi_q_err[kpi])}

    def run_bp_chid(self, nequil=150, nblocks=30, bp=20, ortho=10, pc=10):
        """Unequal-time singlet PAIRING susceptibility (s-wave and d_{x^2-y^2}) for
        the single-band square lattice: C_a(tau_l) = <Delta_a(tau_l) Delta_a^dag(0)>
        for tau_l = l*dt and chi_a = integral C_a(tau) dtau (trapezoidal). Reports
        the FULL susceptibility, the BUBBLE (disconnected, from the ensemble-averaged
        time-displaced GFs: C^0_a(tau)=sum P^up_avg (F_a P^dn_avg F_a^T)), and the
        VERTEX = full - bubble (the interaction-induced connected pairing). Per-block
        error bars. C_a(tau_0) is the equal-time pair structure factor."""
        Fs, Fd, _ = self._bond_factors()
        Ff = {"s": Fs, "d": Fd}
        for it in range(nequil):
            self.step()
            if (it + 1) % ortho == 0: self.reorthogonalize()
            if (it + 1) % pc == 0: self.pop_control()
            self._maybe_update_trial(it)
        # per block, per channel: full and vertex C(tau) arrays
        acc = {f"{t}_{k}": [] for t in ("s", "d") for k in ("full", "vertex")}
        for blk in range(nblocks):
            self.reorthogonalize()
            ket_up = self.walkers.phi_up.copy(); ket_dn = self.walkers.phi_dn.copy()
            rec = [[] for _ in range(self.nw)]
            for _ in range(bp):
                self.prop.step_record(self.walkers, rec)
            Cs, Cd, Pu_s, Pd_s, W = self.est.chid_block(
                self.walkers, ket_up, ket_dn, rec, bp, Fs, Fd)
            if W <= 0:
                self.pop_control(); continue
            Pu = Pu_s / W; Pd = Pd_s / W                 # ensemble-averaged GFs -> bubble
            for tag, Cfull in (("s", Cs / W), ("d", Cd / W)):
                F = Ff[tag]
                bubble = np.array([float(np.sum(Pu[l] * (F @ Pd[l] @ F.T)))
                                   for l in range(bp + 1)])
                acc[f"{tag}_full"].append(Cfull)
                acc[f"{tag}_vertex"].append(Cfull - bubble)
            self.pop_control()
        taus = self.dt * np.arange(bp + 1)
        out = {"nsites": self.n, "bp": bp, "dt": self.dt, "taus": taus.tolist()}
        for key, blocks in acc.items():
            Bm = np.array(blocks)
            chi_blk = self.dt * (Bm[:, 1:-1].sum(axis=1) + 0.5 * (Bm[:, 0] + Bm[:, -1]))
            tag, kind = key.split("_")             # e.g. d_vertex
            suffix = "" if kind == "full" else "_vertex"
            out[f"C{tag}tau{suffix}"] = Bm.mean(axis=0).tolist()
            out[f"chi_{tag}{suffix}"] = float(chi_blk.mean())
            out[f"chi_{tag}{suffix}_err"] = float(chi_blk.std() / np.sqrt(len(chi_blk))
                                                  if len(chi_blk) > 1 else 0.0)
        return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=4); ap.add_argument("--ly", type=int, default=4)
    ap.add_argument("--nup", type=int, default=1); ap.add_argument("--ndn", type=int, default=1)
    ap.add_argument("--U", type=float, default=3.0); ap.add_argument("--t", type=float, default=1.0)
    ap.add_argument("--dt", type=float, default=0.01); ap.add_argument("--nw", type=int, default=400)
    ap.add_argument("--nequil", type=int, default=200); ap.add_argument("--nmeas", type=int, default=400)
    ap.add_argument("-o", "--out", help="write common-schema JSON")
    ap.add_argument("--model", choices=["hubbard", "altermagnet"], default="hubbard")
    ap.add_argument("--t1", type=float, default=-1.0); ap.add_argument("--t2", type=float, default=-1.0)
    ap.add_argument("--t3", type=float, default=-1.0); ap.add_argument("--t4", type=float, default=-1.0)
    ap.add_argument("--uxy", type=float, default=0.0); ap.add_argument("--v", type=float, default=0.0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--bp", type=int, default=0,
                    help="back-propagation length (0 = mixed estimator). >0 removes "
                         "the mixed-estimator bias for uxy/v.")
    ap.add_argument("--corr", action="store_true",
                    help="measure equal-time GF + spin/charge correlations (back-prop) and write JSON")
    ap.add_argument("--trial", choices=["fixed", "adaptive", "sample1", "multidet"], default="fixed",
                    help="trial wavefunction: fixed free-electron GS, or adaptive "
                         "self-consistent natural-orbital trial (Phase 6 WS3)")
    ap.add_argument("--trial-every", type=int, default=10,
                    help="adaptive trial update cadence (steps)")
    ap.add_argument("--trial-mix", type=float, default=0.5,
                    help="adaptive trial damping (0=full replace, 1=no update)")
    ap.add_argument("--trial-k", type=int, default=8,
                    help="multidet trial: number of determinants (sampled walkers)")
    a = ap.parse_args()
    tkw = dict(trial=a.trial, trial_every=a.trial_every, trial_mix=a.trial_mix,
               trial_k=a.trial_k)
    if a.model == "altermagnet":
        import sys, os
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "ed"))
        from altermagnet_ed import build_hopping
        K = build_hopping(a.lx, a.ly, a.t1, a.t2, a.t3, a.t4)   # two-orbital, n=2*lx*ly
        qmc = CPMC(a.lx, a.ly, a.nup, a.ndn, U=a.U, dt=a.dt, nwalkers=a.nw, seed=a.seed,
                   K=K, uxy=a.uxy, v=a.v, **tkw)
        nsites = 2 * a.lx * a.ly
        code = "pyqmc-altermagnet"
    else:
        qmc = CPMC(a.lx, a.ly, a.nup, a.ndn, a.t, a.U, a.dt, a.nw, seed=a.seed, **tkw)
        nsites = a.lx * a.ly
        code = "pyqmc-cpmc"
    if a.corr:
        bp = a.bp if a.bp > 0 else 12
        rec = qmc.run_bp_obs(nequil=a.nequil, nblocks=max(a.nmeas // bp, 12), bp=bp)
        import json
        out = a.out or "/dev/stdout"
        with open(out, "w") as f:
            json.dump(rec, f)
        print(f"wrote correlations ({a.lx}x{a.ly}, nsites={rec['nsites']}) -> {a.out or 'stdout'}")
        raise SystemExit(0)
    if a.bp > 0:
        e, err = qmc.run_bp(nequil=a.nequil, nblocks=max(a.nmeas // a.bp, 10), bp=a.bp)
        est = f"bp={a.bp}"
    else:
        e, err = qmc.run(a.nequil, a.nmeas)
        est = "mixed"
    print(f"CPMC(py,{a.model},{est}) E = {e:.4f} +/- {err:.4f}   (lx={a.lx} ly={a.ly} nup={a.nup} ndn={a.ndn} U={a.U})")
    if a.out:
        import json
        rec = {
            "code": code,
            "model": {"nsites": nsites, "ne": a.nup + a.ndn, "nup": a.nup,
                      "ndn": a.ndn, "U": a.U, "t0": a.t, "lx": a.lx, "ly": a.ly},
            "observables": {
                "energy_total": {"value": e, "error": err},
                "energy_per_site": {"value": e / nsites, "error": err / nsites},
            },
        }
        with open(a.out, "w") as f:
            json.dump(rec, f, indent=2)
        print(f"wrote {a.out}")
