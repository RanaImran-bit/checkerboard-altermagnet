#!/usr/bin/env python3
"""Finite-temperature constrained-path AFQMC for the (altermagnet) Hubbard model
(Phase 2 of docs/PLAN_dqmc_cpqmc.md). The third, sign-controlled finite-T method for the
d-wave pairing cross-validation -- complements the bare finite-T DQMC (sign-problem-prone
at low filling) and the T=0 CPQMC.

Formulation (grand canonical, Z = Tr e^{-beta H}):
  e^{-beta H} = prod_l e^{-dt K'} e^{-dt V_l}; particle-hole symmetric discrete Hirsch HS
  of U -> Ising fields x[:,l]; mu_eff = mu - U/2 folded into K'.  B_l^s(x) = expK_s diag(e^{s lam x}).
Each WALKER builds the path propagator M_l^s = B_l...B_1 slice by slice. With the
non-interacting TRIAL propagator b0_s = expK_s (U=0) and its remaining-slice power
  T_l^s = (b0_s)^{L-l},
the importance function is  I_l = prod_s det(I + T_l^s M_l^s).  Then
  I_0 = prod_s det(I + (b0_s)^L) = Z_T   (trial partition fn),
  I_L = prod_s det(I + M_L^s)     = W(x) (exact path weight),
so  prod_l r_l = I_L/I_0 = W(x)/Z_T  and  <prod r_l>_HS = Z/Z_T  -- EXACT (free projection).
CONSTRAINT: drop a walker whenever an incremental ratio r_l <= 0 (Zhang finite-T CP) ->
tames the sign at the cost of a (trial-dependent) bias; releasing it (carry sign) is exact.

Observables: at the end of the path, the equal-time Green's G_s = (I + M_L^s)^{-1}; thermal
average <O> = sum_paths w_path O[G_path] / sum_paths w_path.

Gates (vs finite-T ED, dqmc.ed_finite_T): U=0 exact; free-projection == DQMC/ED; constrained
within bias of ED with <sign>->1 in the controlled regime.

    python code/ftcpqmc_py/ftcpmc.py --lx 2 --ly 2 --U 4 --mu 2 --beta 3 --ed
"""
from __future__ import annotations
import os, sys, argparse
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "dqmc_py"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "pyqmc"))
from dqmc import bond_factors, reduce_mat, _shift_index, ed_finite_T  # reuse model + ED + reduction
from cpqmc import am_hopping


def _expm(M):
    w, V = np.linalg.eigh(M)
    return (V * np.exp(w)) @ V.conj().T


class FTCPMC:
    def __init__(self, lx, ly, U, mu, beta, dt=0.125, tam=0.0, t1=0.0, seed=1,
                 nw=40, nstab=8, constrained=True, stab=False, Kmat=None, uxy=0.0, v=0.0):
        # Kmat: optional explicit (spin-independent) hopping matrix overriding
        # am_hopping -- pass ed/altermagnet_ed.build_hopping(...) for the two-orbital
        # (d_xz,d_yz) altermagnet (n = 2*lx*ly, both spins see Kmat). On-site U only.
        self.lx, self.ly = lx, ly
        if Kmat is None:
            Ku, Kd = am_hopping(lx, ly, 1.0, tam, t1)
            self.n = lx * ly
        else:
            Ku = Kd = np.asarray(Kmat, float)
            self.n = Ku.shape[0]
        self.norb = self.n // (lx * ly)
        self.U, self.beta, self.dt = U, beta, dt
        self.L = int(round(beta / dt)); self.nstab = nstab
        self.nw = nw; self.constrained = constrained; self.stab = stab
        self.rng = np.random.default_rng(seed)
        from dqmc import _dd_terms
        self.uxy, self.v = uxy, v
        self.dd_terms, self.dd_pot, Kshift = _dd_terms(lx, ly, self.norb, dt, uxy, v)
        Ku0 = np.array(Ku, float); Kd0 = np.array(Kd, float)      # ORIGINAL K (physical kinetic)
        self.K = [Ku0, Kd0]; self._Ku0 = Ku0
        Kus = Ku0.copy(); Kds = Kd0.copy()                        # PROPAGATOR K = orig + dd Hartree shift
        Kus[np.diag_indices(self.n)] += Kshift[0]; Kds[np.diag_indices(self.n)] += Kshift[1]
        mu_eff = mu - 0.5 * U
        self.expK = [np.real(_expm(-dt * (Kus - mu_eff * np.eye(self.n)))),
                     np.real(_expm(-dt * (Kds - mu_eff * np.eye(self.n))))]
        self.expK_inv = [np.real(_expm(dt * (Kus - mu_eff * np.eye(self.n)))),
                         np.real(_expm(dt * (Kds - mu_eff * np.eye(self.n))))]
        self.lam = np.arccosh(np.exp(0.5 * dt * abs(U))) if U != 0 else 0.0
        self.sgn = [+1.0, -1.0]
        self.Fs, self.Fd = bond_factors(lx, ly, self.norb)
        # trial propagator powers T_l^s = (b0_s)^{L-l}, l=0..L (b0_s = expK_s, the U=0 B-matrix)
        self.Tpow = [self._powers(self.expK[s]) for s in (0, 1)]
        if stab:                                   # reuse the brute-validated ASvQRD green() of DQMC
            from dqmc import DQMC                    # same uxy/v -> matching expK + dd_terms (s-form)
            self._dq = DQMC(lx, ly, U, mu, beta, dt, tam, t1, seed, nstab=nstab,
                            Kmat=(None if Kmat is None else Ku0), uxy=uxy, v=v)

    def _green_stab(self, xpath, s):
        """Stabilized equal-time G_s=(I+M_L^s)^{-1} from the sampled path via DQMC.green
        (ASvQRD, brute-validated to 1e-13) -- keeps the measurement accurate at large beta."""
        self._dq.x = np.array(xpath).T            # (n, L): column l = slice-(l+1) fields
        return self._dq.green(s, 0)

    def _powers(self, b0):
        """list P[k] = b0^k for k=0..L (P[L-l] = T_l)."""
        P = [np.eye(self.n)]
        for _ in range(self.L):
            P.append(b0 @ P[-1])
        return P

    def _Bvec(self, s, xcol):
        return np.exp(self.sgn[s] * self.lam * xcol)

    def _imp(self, M):
        """log|I_l| and sign: I_l = prod_s det(I + T_l^s M^s). Returns (logabs, sign)."""
        # caller passes T_l^s @ M^s already? here M is dict s-> (T_l^s @ M^s). Use slogdet(I+.)
        la = 0.0; sg = 1.0
        for s in (0, 1):
            sign, ld = np.linalg.slogdet(np.eye(self.n) + M[s])
            la += ld; sg *= sign
        return la, sg

    def _slice_forcebias(self, l, M):
        """Force-biased (heat-bath) sampling of slice-l HS fields, using the importance
        Green's. Returns (logw, sign, alive). Updates M in place to M_l = expK D M_{l-1}.
        A^s = T_l^s expK^s, C^s = M_{l-1}^s; importance det(I + A^s D^s C^s). Sites sampled
        sequentially from a +1 reference via Sherman-Morrison; walker weight gets the
        per-site normalization (1+|rho|)/2; constraint drops the walker if a chosen flip
        carries rho<0 (the per-site CP boundary)."""
        n = self.n
        A = [self.Tpow[s][self.L - l] @ self.expK[s] for s in (0, 1)]
        C = [M[s] for s in (0, 1)]
        # reference d=1 (no HS) matches I_{l-1} = det(I + A^s C^s); Phi^s = (I + A^s C^s)^{-1}
        Phi = []
        for s in (0, 1):
            try:
                Phi.append(np.linalg.inv(np.eye(n) + A[s] @ C[s]))
            except np.linalg.LinAlgError:
                return 0.0, 1.0, False
        # diagonal HS factor for the two field values, per spin, relative to d=1
        dval = {+1: [np.exp(self.lam), np.exp(-self.lam)],    # x=+1: d^up=e^{+lam}, d^dn=e^{-lam}
                -1: [np.exp(-self.lam), np.exp(self.lam)]}     # x=-1
        x = np.zeros(n); logw = 0.0; sign = 1.0
        for i in range(n):
            facs = {}
            for xv in (+1, -1):
                f = 1.0; upd = []
                for s in (0, 1):
                    dd = dval[xv][s] - 1.0                    # delta vs reference d=1
                    u = A[s][:, i] * dd; v = C[s][i, :]
                    fs = 1.0 + v @ (Phi[s] @ u)
                    f *= fs; upd.append((u, v, fs))
                facs[xv] = (f, upd)
            fp = facs[+1][0]; fm = facs[-1][0]
            if self.constrained:                              # CP: zero out negative-weight branches
                ap = fp if fp > 0 else 0.0; am = fm if fm > 0 else 0.0
            else:
                ap = abs(fp); am = abs(fm)
            tot = ap + am
            if tot <= 0.0:
                return 0.0, 1.0, False                        # both branches killed -> drop walker
            logw += np.log(0.5 * (abs(fp) + abs(fm)))         # importance normalization (1/2 HS prior)
            xv = +1 if self.rng.random() < ap / tot else -1
            f, upd = facs[xv]
            if f < 0:
                sign *= -1.0
            x[i] = xv
            for s in (0, 1):                                  # SM update Phi^s for the chosen field
                u, v, fs = upd[s]
                Phu = Phi[s] @ u; vPh = v @ Phi[s]
                Phi[s] = Phi[s] - np.outer(Phu, vPh) / fs
        for s in (0, 1):                                      # M_l^s = expK^s diag(D^s) M_{l-1}^s
            Dvec = np.exp(self.sgn[s] * self.lam * x)
            M[s] = (self.expK[s] * Dvec[None, :]) @ C[s]
        return logw, sign, True, x

    def run_fb_stab(self, nmeas=400, kres=False, chi=False):
        """STABILIZED force-biased CP-DQMC (P0.1). Per slice l, get the stable equal-time
        Green's G_l (past = sampled B's, future = trial b0's) via DQMC.green (ASvQRD), then
        per-site heat-bath with the LOCAL DQMC ratio R_s = 1 + (1-G_l[i,i]) dv (O(1), stable)
        -- avoids the overflowing (b0)^{L-l} of run_fb. Reference field x=0 (HS factor 1 = b0)
        matches I_{l-1}. Valid at large N*beta where run_fb breaks."""
        n, L = self.n, self.L
        dq = self._dq; ar = np.arange(n)
        accW = accE = accN = signsum = abssum = 0.0
        Mc_d = np.zeros((n, n)); Mc_s = np.zeros((n, n)); Mx_d = np.zeros((n, n)); Mx_s = np.zeros((n, n))
        self._AcU = np.zeros((n, n)); self._AcD = np.zeros((n, n))   # bubble: avg equal-time Green's
        self._accSAM = 0.0                                           # altermagnetic-order SF accumulator
        self._AtU = np.zeros((L, n, n)) if chi else None             # bubble: avg time-displaced Green's
        self._AtD = np.zeros((L, n, n)) if chi else None
        if getattr(self, "paireig", False):                          # k-space pairing-matrix accumulators
            from pair_eig import _dft, _neg_k_index
            self._W, self._kxf, self._kyf = _dft(self.lx, self.ly)
            self._idxm = _neg_k_index(self.lx, self.ly)
            self._Pf = np.zeros((n, n), complex)
            self._Acuk = np.zeros((n, n), complex); self._Acdk = np.zeros((n, n), complex)
        if getattr(self, "paireig_tau", False):                      # tau-integrated k-pairing accumulators
            from pair_eig import _dft, _neg_k_index
            self._W, self._kxf, self._kyf = _dft(self.lx, self.ly)
            self._idxm = _neg_k_index(self.lx, self.ly)
            self._Pft = np.zeros((n, n), complex)
            self._AGuk = [np.zeros((n, n), complex) for _ in range(L)]
            self._AGdk = [np.zeros((n, n), complex) for _ in range(L)]
        for _ in range(nmeas):
            for _w in range(self.nw):
                xpath = []; xddpath = []; logw = 0.0; sign = 1.0; alive = True; G = None
                for l in range(1, L + 1):
                    if G is None or (l - 1) % self.nstab == 0:   # re-stabilize from scratch
                        dq.x[:, :] = 0.0                          # future slices = reference (b0)
                        dq.xdd[:, :] = 0.0                        # future dd = reference (s=0 -> factor 1)
                        for k in range(l - 1):
                            dq.x[:, k] = xpath[k]                 # past slices = sampled
                            if self.dd_terms: dq.xdd[:, k] = xddpath[k]
                        G = [dq.green(s, l - 1) for s in (0, 1)]  # stable Green's at slice l-1
                    # else: G was wrapped forward from the previous slice (Fortran sweep pattern)
                    xl = np.zeros(n)
                    for i in range(n):
                        facs = {}
                        for xv in (1.0, -1.0):
                            R = 1.0
                            for s in (0, 1):
                                dv = np.exp(self.sgn[s] * self.lam * xv) - 1.0
                                R *= 1.0 + (1.0 - G[s][i, i]) * dv
                            facs[xv] = R
                        fp, fm = facs[1.0], facs[-1.0]
                        if self.constrained:
                            ap = fp if fp > 0 else 0.0; am = fm if fm > 0 else 0.0
                        else:
                            ap = abs(fp); am = abs(fm)
                        tot = ap + am
                        if tot <= 0.0:
                            alive = False; break
                        logw += np.log(0.5 * (abs(fp) + abs(fm)))
                        xv = 1.0 if self.rng.random() < ap / tot else -1.0
                        if facs[xv] < 0:
                            sign *= -1.0
                        xl[i] = xv
                        for s in (0, 1):                      # SM update of the stable local Green's
                            dv = np.exp(self.sgn[s] * self.lam * xv) - 1.0
                            Rs = 1.0 + (1.0 - G[s][i, i]) * dv
                            G[s] = G[s] - np.outer((ar == i) - G[s][:, i], G[s][i, :]) * (dv / Rs)
                    if not alive:
                        break
                    # --- off-site density-density (uxy/v) fields at slice l (s in {-1,+1}; ref s=0) ---
                    xl_dd = np.zeros(len(self.dd_terms))
                    for t, (a, sa, b, sb, lam, siga, sigb) in enumerate(self.dd_terms):
                        facs = {}
                        for sv in (1.0, -1.0):
                            dva = np.exp(lam * siga * sv) - 1.0
                            Ra = 1.0 + (1.0 - G[sa][a, a]) * dva
                            dvb = np.exp(lam * sigb * sv) - 1.0
                            Gbb = G[sb][b, b]
                            if sb == sa:
                                Gbb += G[sa][b, a] * G[sa][a, b] * (dva / Ra)
                            facs[sv] = Ra * (1.0 + (1.0 - Gbb) * dvb)
                        fp, fm = facs[1.0], facs[-1.0]
                        if self.constrained:
                            ap = fp if fp > 0 else 0.0; am = fm if fm > 0 else 0.0
                        else:
                            ap = abs(fp); am = abs(fm)
                        tot = ap + am
                        if tot <= 0.0:
                            alive = False; break
                        logw += np.log(0.5 * (abs(fp) + abs(fm)))
                        sv = 1.0 if self.rng.random() < ap / tot else -1.0
                        if facs[sv] < 0:
                            sign *= -1.0
                        xl_dd[t] = sv
                        dva = np.exp(lam * siga * sv) - 1.0      # apply G updates for chosen sv
                        g = G[sa]; Ra = 1.0 + (1.0 - g[a, a]) * dva
                        G[sa] = g - np.outer((ar == a) - g[:, a], g[a, :]) * (dva / Ra)
                        dvb = np.exp(lam * sigb * sv) - 1.0
                        g = G[sb]; Rb = 1.0 + (1.0 - g[b, b]) * dvb
                        G[sb] = g - np.outer((ar == b) - g[:, b], g[b, :]) * (dvb / Rb)
                    if not alive:
                        break
                    xpath.append(xl); xddpath.append(xl_dd)
                    if l < L:                                 # wrap G: slice (l-1) -> l, apply B_{l-1}
                        for s in (0, 1):                       # B=expK diag(hs); hs includes on-site U + dd
                            hs = np.exp(self.sgn[s] * self.lam * xl)
                            for t, (a, sa, b, sb, lam, siga, sigb) in enumerate(self.dd_terms):
                                f = xl_dd[t]
                                if sa == s: hs[a] *= np.exp(lam * siga * f)
                                if sb == s: hs[b] *= np.exp(lam * sigb * f)
                            B = self.expK[s] * hs[None, :]
                            Binv = self.expK_inv[s] * (1.0 / hs)[:, None]
                            G[s] = B @ G[s] @ Binv
                if not alive:
                    continue
                ws = np.exp(logw) * sign
                dq.x[:, :] = np.array(xpath).T                # full sampled path
                if self.dd_terms: dq.xdd[:, :] = np.array(xddpath).T
                Gm = [dq.green(s, 0) for s in (0, 1)]         # stable boundary Green's
                accW, accE, accN, signsum, abssum, Mc_d, Mc_s, Mx_d, Mx_s = self._accum(
                    Gm, ws, sign, xpath, chi, accW, accE, accN, signsum, abssum, Mc_d, Mc_s, Mx_d, Mx_s, xddpath)
        return self._finalize(accW, accE, accN, signsum, abssum, Mc_d, Mc_s, Mx_d, Mx_s, kres, chi)

    def _accum(self, G, ws, sign, xpath, chi, accW, accE, accN, signsum, abssum, Mc_d, Mc_s, Mx_d, Mx_s, xddpath=None):
        """Sign-weighted accumulation of energy/density/pairing from a path's Green's."""
        n = self.n
        cdc_u = np.eye(n) - G[0].T; cdc_d = np.eye(n) - G[1].T
        nu = 1.0 - np.diag(G[0]); nd = 1.0 - np.diag(G[1])
        ek = float(np.sum(self.K[0] * cdc_u.T) + np.sum(self.K[1] * cdc_d.T))
        ev = float(self.U * np.sum(nu * nd))
        for (a, sa, b, sb, V) in self.dd_pot:                   # uxy/v potential (Wick from G)
            na = nu if sa == 0 else nd; nb = nu if sb == 0 else nd
            if sa == sb:
                cdc = cdc_u if sa == 0 else cdc_d
                ev += float(V * (na[a] * nb[b] - cdc[a, b] * cdc[b, a]))
            else:
                ev += float(V * na[a] * nb[b])
        accW += ws; accE += ws * (ek + ev); accN += ws * (nu.sum() + nd.sum())
        signsum += sign; abssum += 1.0
        lxy = self.lx * self.ly; gam = np.where(np.arange(n) < lxy, 1.0, -1.0)   # AM order SF
        sz = 0.25 * (np.outer(nu - nd, nu - nd) + np.diag(nu + nd) - cdc_u * cdc_u.T - cdc_d * cdc_d.T)
        self._accSAM += ws * float(gam @ sz @ gam) / n
        Mc_d += ws * (cdc_u * (self.Fd @ cdc_d @ self.Fd.T))
        Mc_s += ws * (cdc_u * (self.Fs @ cdc_d @ self.Fs.T))
        self._AcU += ws * cdc_u; self._AcD += ws * cdc_d        # bubble: averaged equal-time Green's
        if getattr(self, "paireig", False):                     # k-space pairing matrix P(k,k')=cu(k,k')cd(-k,-k')
            W = self._W; cuk = W @ cdc_u @ W.conj().T; cdk = W @ cdc_d @ W.conj().T
            cdk_m = cdk[np.ix_(self._idxm, self._idxm)]
            self._Pf += ws * (cuk * cdk_m) / n
            self._Acuk += ws * cuk; self._Acdk += ws * cdk_m
        if chi:
            dq = self._dq; dq.x[:, :] = np.array(xpath).T    # path on the DQMC for stable _green_tau
            if self.dd_terms and xddpath is not None: dq.xdd[:, :] = np.array(xddpath).T
            Gt = [G[0].copy(), G[1].copy()]; xd = np.zeros((n, n)); xs = np.zeros((n, n))
            for l in range(self.L):
                if l > 0:
                    if l % dq.nstab == 0:                    # Phase-1 stable time-displaced (avoid raw blow-up)
                        Gt = [dq._green_tau(0, l), dq._green_tau(1, l)]
                    else:
                        for s in (0, 1):
                            Bv = np.exp(self.sgn[s] * self.lam * xpath[l - 1])
                            for t, (a, sa, b, sb, lam, siga, sigb) in enumerate(self.dd_terms):
                                f = xddpath[l - 1][t]
                                if sa == s: Bv[a] *= np.exp(lam * siga * f)
                                if sb == s: Bv[b] *= np.exp(lam * sigb * f)
                            Gt[s] = (self.expK[s] * Bv[None, :]) @ Gt[s]
                gu, gd = Gt
                xd += gu * (self.Fd @ gd @ self.Fd.T); xs += gu * (self.Fs @ gd @ self.Fs.T)
                self._AtU[l] += ws * gu; self._AtD[l] += ws * gd  # bubble: avg time-displaced Green's
                if getattr(self, "paireig_tau", False):           # tau-integrated k-space pairing matrix
                    Guk = self._W @ gu @ self._W.conj().T; Gdk = self._W @ gd @ self._W.conj().T
                    Gdk_m = Gdk[np.ix_(self._idxm, self._idxm)]
                    self._Pft += ws * self.dt * (Guk * Gdk_m)
                    self._AGuk[l] += ws * Guk; self._AGdk[l] += ws * Gdk_m
            Mx_d += ws * self.dt * xd; Mx_s += ws * self.dt * xs
        return accW, accE, accN, signsum, abssum, Mc_d, Mc_s, Mx_d, Mx_s

    def _finalize(self, accW, accE, accN, signsum, abssum, Mc_d, Mc_s, Mx_d, Mx_s, kres, chi):
        out = dict(energy=accE / accW, dens=accN / accW / self.n, sign=signsum / abssum,
                   S_AM=float(self._accSAM / accW),
                   Sd=float(Mc_d.sum() / accW), Ss=float(Mc_s.sum() / accW),
                   chid=float(Mx_d.sum() / accW), chis=float(Mx_s.sum() / accW), npaths=int(abssum))
        if getattr(self, "_AcU", None) is not None:          # scalar equal-time VERTEX (q=0)
            acu = self._AcU / accW; acd = self._AcD / accW
            out["SdV"] = float((Mc_d / accW - acu * (self.Fd @ acd @ self.Fd.T)).sum())
            out["SsV"] = float((Mc_s / accW - acu * (self.Fs @ acd @ self.Fs.T)).sum())
        if chi:                                              # scalar tau-integrated VERTEX (q=0)
            bubd = sum((self._AtU[l] / accW) * (self.Fd @ (self._AtD[l] / accW) @ self.Fd.T)
                       for l in range(self.L))
            bubs = sum((self._AtU[l] / accW) * (self.Fs @ (self._AtD[l] / accW) @ self.Fs.T)
                       for l in range(self.L))
            out["chidV"] = float((Mx_d / accW - self.dt * bubd).sum())
            out["chisV"] = float((Mx_s / accW - self.dt * bubs).sum())
        if getattr(self, "paireig_tau", False):                  # tau-integrated pairing eigenvalue
            Pc = self._Pft / accW - sum(self.dt * (self._AGuk[l] / accW) * (self._AGdk[l] / accW)
                                        for l in range(self.L))
            Pc = 0.5 * (Pc + Pc.conj().T); w, v = np.linalg.eigh(Pc); phi = v[:, -1]
            fd = np.cos(self._kxf) - np.cos(self._kyf); fs = np.cos(self._kxf) + np.cos(self._kyf)
            out["pe_lam"] = float(w[-1]); out["pe_lam2"] = float(w[-2])
            out["pe_dov"] = float(abs(np.vdot(fd / np.linalg.norm(fd), phi)))
            out["pe_sov"] = float(abs(np.vdot(fs / np.linalg.norm(fs), phi)))
        if kres:
            shift = _shift_index(self.lx, self.ly)
            out["corr_d"] = reduce_mat(Mc_d / accW, shift); out["corr_s"] = reduce_mat(Mc_s / accW, shift)
            acu = self._AcU / accW; acd = self._AcD / accW    # VERTEX = FULL - BUBBLE (matches DQMC)
            out["corrV_d"] = reduce_mat(Mc_d / accW - acu * (self.Fd @ acd @ self.Fd.T), shift)
            out["corrV_s"] = reduce_mat(Mc_s / accW - acu * (self.Fs @ acd @ self.Fs.T), shift)
            if chi:
                out["susc_d"] = reduce_mat(Mx_d / accW, shift); out["susc_s"] = reduce_mat(Mx_s / accW, shift)
                bubd = sum((self._AtU[l] / accW) * (self.Fd @ (self._AtD[l] / accW) @ self.Fd.T)
                           for l in range(self.L))
                bubs = sum((self._AtU[l] / accW) * (self.Fs @ (self._AtD[l] / accW) @ self.Fs.T)
                           for l in range(self.L))
                out["suscV_d"] = reduce_mat(Mx_d / accW - self.dt * bubd, shift)
                out["suscV_s"] = reduce_mat(Mx_s / accW - self.dt * bubs, shift)
        if getattr(self, "paireig", False):                   # leading connected pairing eigenvalue
            Pc = self._Pf / accW - (self._Acuk / accW) * (self._Acdk / accW) / self.n
            Pc = 0.5 * (Pc + Pc.conj().T)
            w, v = np.linalg.eigh(Pc); phi = v[:, -1]
            fd = np.cos(self._kxf) - np.cos(self._kyf); fs = np.cos(self._kxf) + np.cos(self._kyf)
            out["lam"] = float(w[-1]); out["lam2"] = float(w[-2])
            out["d_overlap"] = float(abs(np.vdot(fd / np.linalg.norm(fd), phi)))
            out["s_overlap"] = float(abs(np.vdot(fs / np.linalg.norm(fs), phi)))
        return out

    def run_fb(self, nmeas=400, kres=False, chi=False):
        """Force-biased finite-T AFQMC: a fresh importance-sampled path per measurement.
        chi=True also measures the tau-integrated d/s pair susceptibility via the
        time-displaced Green's G_s(tau_l,0)=B(l,0)G_s(0) along the sampled path (= DQMC)."""
        n, L = self.n, self.L
        accW = 0.0; accE = 0.0; accN = 0.0; signsum = 0.0; abssum = 0.0
        Mc_d = np.zeros((n, n)); Mc_s = np.zeros((n, n))
        Mx_d = np.zeros((n, n)); Mx_s = np.zeros((n, n))
        AcU = np.zeros((n, n)); AcD = np.zeros((n, n))            # <c^+c> avg for the equal-time bubble
        AtU = np.zeros((L, n, n)) if chi else None               # <G(l)> avg for the susc bubble
        AtD = np.zeros((L, n, n)) if chi else None
        for _ in range(nmeas):
            for _w in range(self.nw):
                M = [np.eye(n), np.eye(n)]
                logw = 0.0; sign = 1.0; alive = True; xpath = []
                for l in range(1, L + 1):
                    dlw, dsg, ok, xl = self._slice_forcebias(l, M)
                    if not ok:
                        alive = False; break
                    logw += dlw; sign *= dsg; xpath.append(xl)
                if not alive:
                    continue
                ws = np.exp(logw) * sign                       # importance weight x sign (Z_T scale cancels)
                if self.stab:                                  # stabilized measurement Green's
                    G = [self._green_stab(xpath, s) for s in (0, 1)]
                else:
                    G = [np.linalg.solve(np.eye(n) + M[s], np.eye(n)) for s in (0, 1)]
                cdc_u = np.eye(n) - G[0].T; cdc_d = np.eye(n) - G[1].T
                nu = 1.0 - np.diag(G[0]); nd = 1.0 - np.diag(G[1])
                ek = float(np.sum(self.K[0] * cdc_u.T) + np.sum(self.K[1] * cdc_d.T))
                ev = float(self.U * np.sum(nu * nd))
                accW += ws; accE += ws * (ek + ev); accN += ws * (nu.sum() + nd.sum())
                signsum += sign; abssum += 1.0
                Mc_d += ws * (cdc_u * (self.Fd @ cdc_d @ self.Fd.T))
                Mc_s += ws * (cdc_u * (self.Fs @ cdc_d @ self.Fs.T))
                AcU += ws * cdc_u; AcD += ws * cdc_d           # bubble: averaged equal-time Green's
                if chi:                                        # tau-integrated susceptibility
                    Gt = [G[0].copy(), G[1].copy()]            # G_s(tau_0,0)=<c c^dag> equal-time
                    xd = np.zeros((n, n)); xs = np.zeros((n, n))
                    for l in range(L):
                        if l > 0:                              # propagate by B_{l-1}=expK diag(HS_{l-1})
                            for s in (0, 1):
                                Bv = np.exp(self.sgn[s] * self.lam * xpath[l - 1])
                                Gt[s] = (self.expK[s] * Bv[None, :]) @ Gt[s]
                        gu, gd = Gt
                        xd += gu * (self.Fd @ gd @ self.Fd.T); xs += gu * (self.Fs @ gd @ self.Fs.T)
                        AtU[l] += ws * gu; AtD[l] += ws * gd   # bubble: averaged time-displaced Green's
                    Mx_d += ws * self.dt * xd; Mx_s += ws * self.dt * xs
        out = dict(energy=accE / accW, dens=accN / accW / n, sign=signsum / abssum,
                   Sd=float(Mc_d.sum() / accW), Ss=float(Mc_s.sum() / accW),
                   chid=float(Mx_d.sum() / accW), chis=float(Mx_s.sum() / accW), npaths=int(abssum))
        if kres:
            shift = _shift_index(self.lx, self.ly)
            out["corr_d"] = reduce_mat(Mc_d / accW, shift); out["corr_s"] = reduce_mat(Mc_s / accW, shift)
            acu = AcU / accW; acd = AcD / accW                # VERTEX = FULL - BUBBLE (matches DQMC)
            out["corrV_d"] = reduce_mat(Mc_d / accW - acu * (self.Fd @ acd @ self.Fd.T), shift)
            if chi:
                out["susc_d"] = reduce_mat(Mx_d / accW, shift); out["susc_s"] = reduce_mat(Mx_s / accW, shift)
                bub = sum((AtU[l] / accW) * (self.Fd @ (AtD[l] / accW) @ self.Fd.T) for l in range(L))
                out["suscV_d"] = reduce_mat(Mx_d / accW - self.dt * bub, shift)
        return out

    def run(self, nmeas=2000, kres=False):
        n, L, nw = self.n, self.L, self.nw
        # walkers: path propagators M[s] (n x n), weight w, sign carry
        accW = 0.0; accE = 0.0; accN = 0.0; signsum = 0.0; abssum = 0.0
        Mc_d = np.zeros((n, n)); Mc_s = np.zeros((n, n))
        for _ in range(nmeas):
            # one full imaginary-time path per "measurement" (a fresh population of 1 path here;
            # we draw nw independent paths per meas and average)
            for _w in range(nw):
                M = [np.eye(n), np.eye(n)]
                # I_0 = prod_s det(I + (b0)^L M) with M=I -> det(I + Tpow[L])
                Il_log, Il_sg = self._imp([self.Tpow[s][L] for s in (0, 1)])
                wlog = 0.0; wsg = 1.0; alive = True
                for l in range(1, L + 1):
                    x = self.rng.choice([-1.0, 1.0], size=n)
                    for s in (0, 1):
                        M[s] = (self.expK[s] * self._Bvec(s, x)[None, :]) @ M[s]
                    Tl = [self.Tpow[s][L - l] @ M[s] for s in (0, 1)]
                    new_log, new_sg = self._imp(Tl)
                    # incremental ratio r_l = I_l / I_{l-1}
                    r_log = new_log - Il_log; r_sg = new_sg * Il_sg
                    if self.constrained and r_sg <= 0:
                        alive = False; break
                    wlog += r_log; wsg *= r_sg
                    Il_log, Il_sg = new_log, new_sg
                    # NOTE: scale-preserving QR stabilization of M needed for large beta/lattice
                    # (TODO); at small beta the raw product stays in double-precision range.
                if not alive:
                    continue
                w = np.exp(wlog - 0.0) * (wsg if not self.constrained else 1.0)
                wabs = np.exp(wlog); sg = wsg
                # equal-time Green's from the (rescaled) full path: G_s = (I + M_L^s)^{-1}
                G = [np.linalg.solve(np.eye(n) + M[s], np.eye(n)) for s in (0, 1)]
                cdc_u = np.eye(n) - G[0].T; cdc_d = np.eye(n) - G[1].T
                nu = 1.0 - np.diag(G[0]); nd = 1.0 - np.diag(G[1])
                ek = float(np.sum(self.K[0] * cdc_u.T) + np.sum(self.K[1] * cdc_d.T))
                ev = float(self.U * np.sum(nu * nd))
                wsigned = wabs * sg
                accW += wsigned; accE += wsigned * (ek + ev); accN += wsigned * (nu.sum() + nd.sum())
                signsum += sg; abssum += 1.0
                Mc_d += wsigned * (cdc_u * (self.Fd @ cdc_d @ self.Fd.T))
                Mc_s += wsigned * (cdc_u * (self.Fs @ cdc_d @ self.Fs.T))
        out = dict(energy=accE / accW, dens=accN / accW / n, sign=signsum / abssum,
                   Sd=float(Mc_d.sum() / accW), Ss=float(Mc_s.sum() / accW), npaths=int(abssum))
        if kres:
            shift = _shift_index(self.lx, self.ly)
            out["corr_d"] = reduce_mat(Mc_d / accW, shift)
            out["corr_s"] = reduce_mat(Mc_s / accW, shift)
        return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=2); ap.add_argument("--ly", type=int, default=2)
    ap.add_argument("--U", type=float, default=4.0); ap.add_argument("--mu", type=float, default=2.0)
    ap.add_argument("--beta", type=float, default=3.0); ap.add_argument("--dt", type=float, default=0.125)
    ap.add_argument("--tam", type=float, default=0.0); ap.add_argument("--t1", type=float, default=0.0)
    ap.add_argument("--nw", type=int, default=40); ap.add_argument("--nmeas", type=int, default=400)
    ap.add_argument("--seed", type=int, default=1); ap.add_argument("--ed", action="store_true")
    ap.add_argument("--free", action="store_true", help="release constraint (exact, sign-carrying)")
    ap.add_argument("--fb", action="store_true", help="force-biased (heat-bath) importance sampling")
    ap.add_argument("--chi", action="store_true", help="also measure tau-integrated pair susceptibility (fb only)")
    ap.add_argument("--kres", action="store_true")
    ap.add_argument("--stab", action="store_true", help="stabilized (ASvQRD) measurement Green's")
    a = ap.parse_args()
    q = FTCPMC(a.lx, a.ly, a.U, a.mu, a.beta, a.dt, a.tam, a.t1, a.seed,
               nw=a.nw, constrained=not a.free, stab=a.stab)
    if a.fb and a.stab:
        r = q.run_fb_stab(a.nmeas, kres=a.kres, chi=a.chi)
    elif a.fb:
        r = q.run_fb(a.nmeas, kres=a.kres, chi=a.chi)
    else:
        r = q.run(a.nmeas, kres=a.kres)
    mode = "FREE" if a.free else "CONSTRAINED"
    print(f"# FT-CPMC {a.lx}x{a.ly} U={a.U} mu={a.mu} beta={a.beta} (L={q.L}) tam={a.tam} t1={a.t1} [{mode}]")
    print(f"  <sign> = {r['sign']:.4f}   npaths={r['npaths']}")
    print(f"  density = {r['dens']:.5f}")
    print(f"  energy(K+U) = {r['energy']:.5f}")
    print(f"  S_d (q=0) = {r['Sd']:.4f}   S_s = {r['Ss']:.4f}")
    if a.chi and "chid" in r:
        print(f"  chi_d (q=0) = {r['chid']:.4f}   chi_s = {r['chis']:.4f}")
    if a.kres and "corrV_d" in r:
        cv = r["corrV_d"]; print(f"  corrV_d (vtx) maxk={cv['maxk']:.4f} k0={cv['k0']:.4f}")
        if a.chi and "suscV_d" in r:
            xv = r["suscV_d"]; print(f"  suscV_d (vtx) maxk={xv['maxk']:.4f} k0={xv['k0']:.4f}")
    if a.ed:
        e, d, Sd, Ss, chid, chis = ed_finite_T(a.lx, a.ly, a.U, a.mu, a.beta, a.tam, a.t1)
        print(f"  ED  density = {d:.5f}   energy = {e:.5f}   S_d = {Sd:.4f}   S_s = {Ss:.4f}")
        print(f"  ED  chi_d = {chid:.4f}   chi_s = {chis:.4f}")


if __name__ == "__main__":
    main()
