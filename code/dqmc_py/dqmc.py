#!/usr/bin/env python3
"""Focused numpy finite-temperature determinant QMC (BSS) for the (altermagnet) Hubbard
model -- the finite-T benchmark for the CPQMC d-wave pairing susceptibility (Phase 1 of
docs/PLAN_dqmc_cpqmc.md). Kept deliberately small + heavily ED-gated.

Model: H = sum_sigma c^dag K_sigma c + U sum_i n_iu n_id - mu sum_i n_i, with the
spin-dependent altermagnet hopping K_sigma = am_hopping(.,tam,t1) (K = -t convention).

BSS: beta = NT*dt; Trotter e^{-beta H} = prod_l e^{-dt K'} e^{-dt V_l}; particle-hole
symmetric discrete Hirsch HS of U (n_u-1/2)(n_d-1/2) -> Ising fields x[i,l] in {+-1};
mu_eff = mu - U/2 folded into K'. Per-slice spin propagator
    B_l^sigma = e^{-dt(K - mu_eff)} diag(e^{sigma * lam * x[:,l]}),  cosh(lam)=e^{dt U/2}.
Equal-time Green's G^sigma = (I + B_NT..B_1)^{-1} (QR-stabilized). Metropolis sweep over
all (i,l) with Sherman-Morrison updates; weight det(I+A^u)det(I+A^d) (real, sign tracked).
Measures <n>, energy, equal-time d-wave/s-wave pair structure factor; <sign>.

    python code/dqmc_py/dqmc.py --lx 2 --ly 2 --U 4 --mu 2 --beta 4 --ed   # gate vs finite-T ED
"""
from __future__ import annotations
import os, sys, argparse
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "pyqmc"))
from cpqmc import am_hopping


def _shift_index(lx, ly):
    """shift[dx,dy] = site map m -> (x+dx, y+dy) (PBC). Site i = x*ly + y. Matches unified_scan."""
    n = lx * ly; xm = np.arange(n) // ly; ym = np.arange(n) % ly
    S = np.empty((lx, ly, n), dtype=int)
    for dx in range(lx):
        for dy in range(ly):
            S[dx, dy] = ((xm + dx) % lx) * ly + (ym + dy) % ly
    return S


def reduce_mat(mat, shift, rmin=2.0):
    """M[m,n] -> S(R)=sum_m M[m,m+R] (lx,ly grid), P(q)=FFT2(S). dict(maxk,k0,r0,rgt).
    IDENTICAL reduction to pyqmc/unified_scan.reduce_mat (so DQMC q-grid == CPQMC q-grid)."""
    lx, ly, n = shape = shift.shape
    rows = np.arange(n)
    Sg = np.array([[mat[rows, shift[dx, dy]].sum() for dy in range(ly)] for dx in range(lx)])
    Pq = np.real(np.fft.fft2(Sg))
    dxg = np.minimum(np.arange(lx), lx - np.arange(lx))[:, None]
    dyg = np.minimum(np.arange(ly), ly - np.arange(ly))[None, :]
    Rmag = np.sqrt(dxg ** 2 + dyg ** 2) * np.ones((lx, ly))
    rgt = (Sg[Rmag > rmin] / n).mean() if np.any(Rmag > rmin) else float("nan")
    return dict(maxk=float(Pq.max()), k0=float(Sg.sum()), r0=float(Sg[0, 0] / n), rgt=float(rgt),
                Pq=Pq, Sg=Sg)


def bond_factors(lx, ly, norb=1):
    """Intra-orbital NN d-wave (+x,-y) and s-wave (+x,+y) bond form-factor matrices.
    For norb>1 (two-orbital altermagnet) the bonds are built within each orbital
    block, index = orb*lxy + x*ly + y (matches ed/altermagnet_ed.build_hopping and
    cpqmc._bond_factors)."""
    lxy = lx * ly; n = norb * lxy
    Fs = np.zeros((n, n)); Fd = np.zeros((n, n))
    def idx(x, y, orb): return orb * lxy + (x % lx) * ly + (y % ly)
    for orb in range(norb):
        for x in range(lx):
            for y in range(ly):
                m = idx(x, y, orb)
                for (dx, dy, fd) in ((1, 0, 1.0), (-1, 0, 1.0), (0, 1, -1.0), (0, -1, -1.0)):
                    j = idx(x + dx, y + dy, orb); Fs[m, j] += 1.0; Fd[m, j] += fd
    return Fs, Fd


def _dd_terms(lx, ly, norb, dt, uxy, v):
    """Off-site density-density interaction term list for the two-orbital model, in the
    DIFFERENCE/SUM channel with an s=0 field-OFF reference (so both DQMC-Metropolis and the
    CP-DQMC trial-based stabilizer can use it). For each mode pair (alpha,beta):
      n_a n_b = 1/2(n_a+n_b) - 1/2(n_a-n_b)^2   -> exact discrete HS:
      V>0 (repulsive): e^{-dtV n_a n_b} = e^{-dtV(n_a+n_b)/2} * 1/2 sum_s e^{lam s (n_a-n_b)},
                       cosh(lam)=e^{dtV/2}; field on a: e^{lam s}, on b: e^{-lam s}.
      V<0 (attractive): e^{-dtV n_a n_b} = e^{+dtV(n_a+n_b)/2} * 1/2 sum_s e^{lam s (n_a+n_b)},
                       cosh(lam)=e^{-dtV/2}; field on a,b: e^{lam s} (same sign).
    The constant e^{-+dtV(n_a+n_b)/2} is folded into K (returned Kshift[spin,site]). Field
    index k in {0,1} -> s=2k-1 in {-1,+1}; s=0 gives factor 1 (the trial reference). Entry
    (a,sa,b,sb,e1,e2,cf) with cf=[1,1]. Inter-orbital U' (uxy) + neighbour V (+v intra/-v inter,
    +x/+y bonds). Empty for the single band (norb=1)."""
    n = norb * lx * ly
    Kshift = np.zeros((2, n)); terms = []; pot = []
    if norb < 2 or (uxy == 0 and v == 0):
        return terms, pot, Kshift
    lxy = lx * ly
    def idx0(ix, iy, orb): return orb * lxy + (ix % lx) * ly + (iy % ly)
    def add(a, sa, b, sb, V):
        # field s in {-1,0,+1}: factor on (a,sa) = e^{lam*siga*s}, on (b,sb) = e^{lam*sigb*s};
        # s=0 -> factor 1 (trial-off reference). diff channel (V>0): sigb=-siga; sum (V<0): sigb=+siga.
        if V > 0:
            lam = np.arccosh(np.exp(0.5 * dt * V)); siga, sigb = 1.0, -1.0
            Kshift[sa][a] += 0.5 * V; Kshift[sb][b] += 0.5 * V
        else:
            lam = np.arccosh(np.exp(-0.5 * dt * V)); siga, sigb = 1.0, 1.0
            Kshift[sa][a] += -0.5 * V; Kshift[sb][b] += -0.5 * V
        terms.append((a, sa, b, sb, lam, siga, sigb)); pot.append((a, sa, b, sb, V))
    if uxy != 0:
        for r in range(lxy):
            for sa in (0, 1):
                for sb in (0, 1):
                    add(r, sa, r + lxy, sb, uxy)
    if v != 0:
        for ix in range(lx):
            for iy in range(ly):
                for (jx, jy) in (((ix + 1) % lx, iy), (ix, (iy + 1) % ly)):
                    if (jx, jy) == (ix, iy):
                        continue
                    for oi in (0, 1):
                        for oj in (0, 1):
                            a = idx0(ix, iy, oi); b = idx0(jx, jy, oj)
                            Vc = v if oi == oj else -v
                            for sa in (0, 1):
                                for sb in (0, 1):
                                    add(a, sa, b, sb, Vc)
    return terms, pot, Kshift


class DQMC:
    def __init__(self, lx, ly, U, mu, beta, dt=0.125, tam=0.0, t1=0.0, seed=1, nstab=8, tp=0.0,
                 Kmat=None, uxy=0.0, v=0.0):
        # Kmat: optional explicit (spin-independent) hopping matrix. If given it
        # OVERRIDES am_hopping -- pass ed/altermagnet_ed.build_hopping(...) to run
        # the two-orbital (d_xz,d_yz) altermagnet (n = 2*lx*ly). Both spins see Kmat.
        self.lx, self.ly = lx, ly
        if Kmat is None:
            Ku, Kd = am_hopping(lx, ly, 1.0, tam, t1, tp)
            self.n = lx * ly
        else:
            Ku = Kd = np.asarray(Kmat, float)
            self.n = Ku.shape[0]
        self.norb = self.n // (lx * ly)
        self.U, self.beta, self.dt = U, beta, dt
        self.NT = int(round(beta / dt)); self.nstab = nstab
        self.rng = np.random.default_rng(seed)
        self.uxy, self.v = uxy, v
        self.dd_terms, self.dd_pot, Kshift = _dd_terms(lx, ly, self.norb, dt, uxy, v)
        Ku0 = np.array(Ku, float); Kd0 = np.array(Kd, float)   # ORIGINAL K (physical kinetic)
        self.K = [Ku0, Kd0]
        Kus = Ku0.copy(); Kds = Kd0.copy()                     # PROPAGATOR K = orig + dd Hartree shift
        Kus[np.diag_indices(self.n)] += Kshift[0]; Kds[np.diag_indices(self.n)] += Kshift[1]
        mu_eff = mu - 0.5 * U
        self.expK = [np.real(_expm(-dt * (Kus - mu_eff * np.eye(self.n)))),
                     np.real(_expm(-dt * (Kds - mu_eff * np.eye(self.n))))]
        self.lam = np.arccosh(np.exp(0.5 * dt * abs(U))) if U != 0 else 0.0
        self.sgn = [+1.0, -1.0]                     # spin sign in e^{sigma lam x}
        self.x = self.rng.choice([-1.0, 1.0], size=(self.n, self.NT))
        self.Fs, self.Fd = bond_factors(lx, ly, self.norb)
        # off-site density-density (uxy, v) HS fields (difference/sum channel), values +-1;
        # on-site U stays in self.x. (s=0 is the trial-off reference, used by CP-DQMC.)
        self.xdd = (self.rng.choice([-1.0, 1.0], size=(len(self.dd_terms), self.NT))
                    if self.dd_terms else np.zeros((0, self.NT)))

    def _Bvec(self, s, l):
        """diagonal HS factor vector for spin s, slice l: on-site U (Hirsch spin channel)
        times all off-site density-density (uxy/v) factors touching spin s at slice l."""
        d = np.exp(self.sgn[s] * self.lam * self.x[:, l])
        for t, (a, sa, b, sb, lam, siga, sigb) in enumerate(self.dd_terms):
            f = self.xdd[t, l]
            if sa == s: d[a] *= np.exp(lam * siga * f)
            if sb == s: d[b] *= np.exp(lam * sigb * f)
        return d

    def _Aprod(self, s):
        """A = B_{NT}..B_1 for spin s (used for from-scratch G)."""
        A = np.eye(self.n)
        for l in range(self.NT):
            A = self.expK[s] @ (self._Bvec(s, l)[:, None] * A)
        return A

    def green(self, s, l0=0):
        """Equal-time Green's at slice l0: G_{l0} = (I + B_{l0-1}..B_0 B_{NT-1}..B_{l0})^{-1}
        via stable ASvQRD (chain = U diag(D) T, big/small-scale split for the inverse).
        l0=0 is the standard G. Used for restabilization mid-sweep."""
        n = self.n
        Uu = np.eye(n); D = np.ones(n); T = np.eye(n)
        for k in range(self.NT):
            l = (l0 + k) % self.NT
            M = (self.expK[s] @ (self._Bvec(s, l)[:, None] * Uu)) * D[None, :]
            Q, R = np.linalg.qr(M)
            Dg = np.diag(R).copy()
            T = (R / Dg[:, None]) @ T          # unit-diagonal upper-tri times accumulated T
            D = Dg; Uu = Q                      # signed scales
        Db = np.where(np.abs(D) > 1.0, D, 1.0)  # big part
        Ds = np.where(np.abs(D) > 1.0, 1.0, D)  # small part
        inner = (Uu.T / Db[:, None]) + (Ds[:, None] * T)
        return np.linalg.solve(inner, Uu.T / Db[:, None])

    def _udv_chain(self, s, a, b):
        """B(b,a) = B_{b-1}..B_a = U diag(d) V (QR-stabilized; U orth, V unit-upper-tri).
        b>a; B_k = expK diag(HS_k)."""
        n = self.n; U = np.eye(n); d = np.ones(n); V = np.eye(n)
        for k in range(a, b):
            M = (self.expK[s] @ (self._Bvec(s, k)[:, None] * U)) * d[None, :]
            Q, R = np.linalg.qr(M)
            dg = np.diag(R).copy()
            V = (R / dg[:, None]) @ V; d = dg; U = Q
        return U, d, V

    def _green_tau(self, s, l):
        """STABLE time-displaced Green's G(tau_l,0) = <c(tau_l) c^dag(0)> = B(l,0)(I+B(NT,0))^-1.
        = (B(l,0)^{-1} + B(NT,l))^{-1} via big/small split (docs/STABILIZATION.md):
        INNER = D1b^-1 (U1^T V2^-1) D2b^-1 + D1s (V1 U2) D2s ; G = V2^-1 D2b^-1 INNER^-1 D1s V1."""
        if l == 0:
            return self.green(s, 0)
        U1, D1, V1 = self._udv_chain(s, 0, l)          # B(l,0)
        U2, D2, V2 = self._udv_chain(s, l, self.NT)    # B(NT,l)
        D1b = np.where(np.abs(D1) > 1.0, D1, 1.0); D1s = np.where(np.abs(D1) > 1.0, 1.0, D1)
        D2b = np.where(np.abs(D2) > 1.0, D2, 1.0); D2s = np.where(np.abs(D2) > 1.0, 1.0, D2)
        V2i = np.linalg.inv(V2)
        INNER = (U1.T @ V2i) / D1b[:, None] / D2b[None, :] + D1s[:, None] * (V1 @ U2) * D2s[None, :]
        X = np.linalg.solve(INNER, D1s[:, None] * V1)  # INNER^-1 D1s V1
        return V2i @ (X / D2b[:, None])                # V2^-1 D2b^-1 X

    def sweep(self):
        """One full space-time sweep; returns the configuration sign (+-1)."""
        G = [self.green(0, 0), self.green(1, 0)]
        sign = 1.0
        for l in range(self.NT):
            if l > 0 and l % self.nstab == 0:           # restabilize with slice-l Green's
                G = [self.green(0, l), self.green(1, l)]  # (green(s,l) brute-validated to 1e-13)
            for i in range(self.n):
                # propose x[i,l] -> -x[i,l]; delta in the HS factor for each spin at site i
                R = 1.0
                deltas = []
                for s in (0, 1):
                    dv = np.exp(self.sgn[s] * self.lam * (-2.0 * self.x[i, l])) - 1.0
                    Rs = 1.0 + (1.0 - G[s][i, i]) * dv
                    R *= Rs; deltas.append((dv, Rs))
                if self.rng.random() < min(1.0, abs(R)):
                    self.x[i, l] *= -1.0
                    if R < 0:
                        sign *= -1.0
                    for s in (0, 1):
                        dv, Rs = deltas[s]
                        g = G[s]
                        # Sherman-Morrison: G' = G - (dv/R)(e_i - G[:,i]) G[i,:]
                        G[s] = g - np.outer((np.arange(self.n) == i) - g[:, i], g[i, :]) * (dv / Rs)
            # off-site density-density (uxy/v) HS fields: one per term at this slice.
            # Each flip changes two diagonal factors (site a spin sa, site b spin sb; a!=b);
            # joint det ratio = Ra * Rb (Rb on the a-updated Green's if sb==sa) * bosonic cf.
            ar = np.arange(self.n)
            for t, (a, sa, b, sb, lam, siga, sigb) in enumerate(self.dd_terms):
                f = self.xdd[t, l]                        # propose flip +-1 -> -+1
                dva = np.exp(lam * siga * (-2.0 * f)) - 1.0   # new/old factor - 1 on (a,sa)
                dvb = np.exp(lam * sigb * (-2.0 * f)) - 1.0   # on (b,sb)
                Ra = 1.0 + (1.0 - G[sa][a, a]) * dva
                Gbb = G[sb][b, b]
                if sb == sa:                              # a-update shifts G[sb][b,b]
                    Gbb += G[sa][b, a] * G[sa][a, b] * (dva / Ra)
                Rb = 1.0 + (1.0 - Gbb) * dvb
                R = Ra * Rb
                if self.rng.random() < min(1.0, abs(R)):
                    self.xdd[t, l] = -f
                    if R < 0:
                        sign *= -1.0
                    g = G[sa]                             # rank-1 at (a, sa)
                    G[sa] = g - np.outer((ar == a) - g[:, a], g[a, :]) * (dva / Ra)
                    g = G[sb]                             # rank-1 at (b, sb) on current G[sb]
                    Rb2 = 1.0 + (1.0 - g[b, b]) * dvb
                    G[sb] = g - np.outer((ar == b) - g[:, b], g[b, :]) * (dvb / Rb2)
            # propagate G to next slice: G(l+1) = B_l G(l) B_l^{-1} (B_l = expK diag(HS_l))
            for s in (0, 1):
                B = self.expK[s] * self._Bvec(s, l)[None, :]
                G[s] = B @ G[s] @ np.linalg.inv(B)
        return sign

    def measure(self, G):
        """Equal-time observables from G^sigma. Returns dict."""
        gu, gd = G
        nu = 1.0 - np.diag(gu); nd = 1.0 - np.diag(gd)        # <n_i,sigma>
        dens = float((nu.sum() + nd.sum()) / self.n)
        # <c^+_i c_j>_s = delta_ij - g_s[j,i]
        cdc_u = np.eye(self.n) - gu.T; cdc_d = np.eye(self.n) - gd.T
        ek = float(np.sum(self.K[0] * cdc_u.T) + np.sum(self.K[1] * cdc_d.T))
        ev = float(self.U * np.sum(nu * nd))
        for (a, sa, b, sb, V) in self.dd_pot:            # uxy/v potential (Wick from G)
            na = nu if sa == 0 else nd; nb = nu if sb == 0 else nd
            if sa == sb:                                 # same spin: Pauli exchange reduces it
                cdc = cdc_u if sa == 0 else cdc_d
                ev += float(V * (na[a] * nb[b] - cdc[a, b] * cdc[b, a]))
            else:
                ev += float(V * na[a] * nb[b])
        # equal-time pair correlation MATRIX M_a[m,n] = <Delta_a^dag(m)Delta_a(n)>
        # = cu[m,n] * (F_a cd F_a^T)[m,n]  (Wick, q=0 sum recovers S_a); k-resolve via reduce_mat
        Mcorr_d = cdc_u * (self.Fd @ cdc_d @ self.Fd.T)
        Mcorr_s = cdc_u * (self.Fs @ cdc_d @ self.Fs.T)
        # altermagnetic-order structure factor S_AM = <(sum_i g_i S^z_i)^2>/N, g=+1/-1 by orbital
        lxy = self.lx * self.ly; g = np.where(np.arange(self.n) < lxy, 1.0, -1.0)
        sz = 0.25 * (np.outer(nu - nd, nu - nd) + np.diag(nu + nd) - cdc_u * cdc_u.T - cdc_d * cdc_d.T)
        S_AM = float(g @ sz @ g) / self.n
        return dict(dens=dens, energy=ek + ev, ek=ek, ev=ev, nu=nu, nd=nd, S_AM=S_AM,
                    cdc_u=cdc_u, cdc_d=cdc_d, Mcorr_d=Mcorr_d, Mcorr_s=Mcorr_s)

    def chi_pair(self, s=1.0, AtU=None, AtD=None):
        """tau-integrated q=0 pair susceptibility chi_a = int_0^beta dtau <Delta_a(tau) Delta_a^dag(0)>.
        Wick: <Delta_a(tau)Delta_a^dag(0)> = sum_{mn} Gu(tau)_{mn} [F_a Gd(tau) F_a^T]_{mn},
        time-displaced G_s(tau_l,0) = <c(tau_l) c^dag(0)> = B(l,0) G_s(0). If AtU/AtD given,
        also accumulate s*Gu(l), s*Gd(l) per slice (for the susceptibility BUBBLE = vertex)."""
        Gt = [self.green(0, 0), self.green(1, 0)]   # G_s(tau_0,0) = equal-time <c c^dag>
        Mchi_d = np.zeros((self.n, self.n)); Mchi_s = np.zeros((self.n, self.n))
        for l in range(self.NT):
            if l > 0:                                # propagate tau_{l-1} -> tau_l: B_{l-1} = expK diag(HS_{l-1})
                if l % self.nstab == 0:              # re-stabilize (stable two-sided UDV) -- avoids
                    Gt = [self._green_tau(0, l), self._green_tau(1, l)]   # the raw-product blow-up at large NT
                else:
                    for sp in (0, 1):
                        Gt[sp] = (self.expK[sp] * self._Bvec(sp, l - 1)[None, :]) @ Gt[sp]
            gu, gd = Gt
            Mchi_d += gu * (self.Fd @ gd @ self.Fd.T)   # tau-integrand MATRIX (k-resolvable)
            Mchi_s += gu * (self.Fs @ gd @ self.Fs.T)
            if AtU is not None:
                AtU[l] += s * gu; AtD[l] += s * gd       # bubble averages <Gu(l)>, <Gd(l)>
        return self.dt * Mchi_d, self.dt * Mchi_s       # q=0 sum recovers scalar chi_a

    def run(self, nwarm=200, nmeas=400, nsub=2, chi=False, kres=False):
        for _ in range(nwarm):
            self.sweep()
        n, L = self.n, self.NT
        acc = {"dens": 0.0, "energy": 0.0, "S_AM": 0.0}; sw = 0.0; saw = 0.0
        Mc_d = np.zeros((n, n)); Mc_s = np.zeros((n, n))   # equal-time corr matrices (sign-wtd)
        Mx_d = np.zeros((n, n)); Mx_s = np.zeros((n, n))   # tau-integrated susc matrices
        AcU = np.zeros((n, n)); AcD = np.zeros((n, n))     # <cdc_u>, <cdc_d> for the equal-time bubble
        AtU = np.zeros((L, n, n)) if chi else None         # <Gu(l)>, <Gd(l)> for the susc bubble
        AtD = np.zeros((L, n, n)) if chi else None
        for _ in range(nmeas):
            for _ in range(nsub):
                s = self.sweep()
            G = [self.green(0), self.green(1)]
            m = self.measure(G)
            Mc_d += s * m["Mcorr_d"]; Mc_s += s * m["Mcorr_s"]
            AcU += s * m["cdc_u"]; AcD += s * m["cdc_d"]
            if chi:
                xd, xs = self.chi_pair(s, AtU, AtD)
                Mx_d += s * xd; Mx_s += s * xs
            for k in acc:
                acc[k] += s * m[k]
            sw += s; saw += abs(s)
        Mc_d /= sw; Mc_s /= sw; Mx_d /= sw; Mx_s /= sw; AcU /= sw; AcD /= sw
        # VERTEX = FULL - BUBBLE (bubble = pair kernel of the config-averaged Green's)
        Vc_d = Mc_d - AcU * (self.Fd @ AcD @ self.Fd.T)    # equal-time d-wave vertex
        Vc_s = Mc_s - AcU * (self.Fs @ AcD @ self.Fs.T)    # equal-time ext-s vertex
        out = dict(dens=acc["dens"] / sw, energy=acc["energy"] / sw, S_AM=acc["S_AM"] / sw,
                   Sd=float(Mc_d.sum()), Ss=float(Mc_s.sum()),
                   SdV=float(Vc_d.sum()), SsV=float(Vc_s.sum()),
                   chid=float(Mx_d.sum()), chis=float(Mx_s.sum()),
                   sign=sw / saw, nmeas=nmeas)
        if chi:
            AtU /= sw; AtD /= sw
            Vx_d = Mx_d - self.dt * sum(AtU[l] * (self.Fd @ AtD[l] @ self.Fd.T) for l in range(L))
            Vx_s = Mx_s - self.dt * sum(AtU[l] * (self.Fs @ AtD[l] @ self.Fs.T) for l in range(L))
            out["chidV"] = float(Vx_d.sum()); out["chisV"] = float(Vx_s.sum())
        if kres:                                           # k-resolve via the unified_scan reduction
            shift = _shift_index(self.lx, self.ly)
            out["corr_d"] = reduce_mat(Mc_d, shift); out["corr_s"] = reduce_mat(Mc_s, shift)
            out["corrV_d"] = reduce_mat(Vc_d, shift); out["corrV_s"] = reduce_mat(Vc_s, shift)
            if chi:
                out["susc_d"] = reduce_mat(Mx_d, shift); out["susc_s"] = reduce_mat(Mx_s, shift)
                out["suscV_d"] = reduce_mat(Vx_d, shift); out["suscV_s"] = reduce_mat(Vx_s, shift)
        return out


def _expm(M):
    w, V = np.linalg.eigh(M)
    return (V * np.exp(w)) @ V.conj().T


def ed_finite_T(lx, ly, U, mu, beta, tam=0.0, t1=0.0):
    """Finite-T ED: Tr(e^{-beta H} O)/Tr(e^{-beta H}) for energy and density, full Fock."""
    from quspin.basis import spinful_fermion_basis_general
    from quspin.operators import hamiltonian
    n = lx * ly; basis = spinful_fermion_basis_general(n)
    nc = dict(check_pcon=False, check_symm=False, check_herm=False)
    Ku, Kd = am_hopping(lx, ly, 1.0, tam, t1)
    up = [[Ku[i, j], i, j] for i in range(n) for j in range(n) if abs(Ku[i, j]) > 1e-15]
    dn = [[Kd[i, j], i, j] for i in range(n) for j in range(n) if abs(Kd[i, j]) > 1e-15]
    inter = [[U, i, i] for i in range(n)]
    mut = [[-mu, i] for i in range(n)]
    Hk = hamiltonian([["+-|", up], ["|+-", dn]], [], basis=basis, dtype=np.float64, **nc).toarray()
    Hu = hamiltonian([["n|n", inter]], [], basis=basis, dtype=np.float64, **nc).toarray()
    Hmu = hamiltonian([["n|", mut], ["|n", mut]], [], basis=basis, dtype=np.float64, **nc).toarray()
    H = Hk + Hu + Hmu
    Nop = hamiltonian([["n|", [[1.0, i] for i in range(n)]], ["|n", [[1.0, i] for i in range(n)]]],
                      [], basis=basis, dtype=np.float64, **nc).toarray()
    # equal-time q=0 pair structure factors S_a = Tr(rho Delta_a^dag Delta_a)/Z
    from dqmc import bond_factors
    Fs, Fd = bond_factors(lx, ly)
    def pair_op(F):
        terms = [[F[m, j], m, j] for m in range(n) for j in range(n) if abs(F[m, j]) > 1e-15]
        return hamiltonian([["+|+", terms]], [], basis=basis, dtype=np.float64, **nc).toarray()
    Dd = pair_op(Fd); Ds = pair_op(Fs)
    w, V = np.linalg.eigh(H)
    bw = np.exp(-beta * (w - w.min())); Z = bw.sum()
    Hku = Hk + Hu                         # energy = kinetic + U (no mu term), to match DQMC
    def thermal(O): return float(np.sum(bw * np.real(np.einsum("ik,ij,jk->k", V.conj(), O, V))) / Z)
    e = thermal(Hku); dens = thermal(Nop) / n
    Sd = thermal(Dd @ Dd.conj().T); Ss = thermal(Ds @ Ds.conj().T)
    # static (tau-integrated) pair susceptibility chi_a = int_0^beta dtau <Delta_a(tau)Delta_a^dag(0)>
    # Kubo/Lehmann: chi_a = (1/Z) sum_{ab} |<b|Delta_a^dag|a>|^2 K(E_a,E_b),
    #   K = (e^{-beta E_b}-e^{-beta E_a})/(E_a-E_b),  -> beta e^{-beta E_a} when E_a=E_b.
    def chi(D):
        Mt = V.conj().T @ D @ V                      # Mt[b,a] = <b|Delta^dag|a>
        diff = w[:, None] - w[None, :]               # E_a - E_b  (a=row, b=col)
        small = np.abs(diff) < 1e-9
        K = np.where(small, beta * bw[:, None],
                     (bw[None, :] - bw[:, None]) / np.where(small, 1.0, diff))
        return float((np.abs(Mt).T ** 2 * K).sum() / Z)
    chid = chi(Dd); chis = chi(Ds)
    return e, dens, Sd, Ss, chid, chis


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=2); ap.add_argument("--ly", type=int, default=2)
    ap.add_argument("--U", type=float, default=4.0); ap.add_argument("--mu", type=float, default=2.0)
    ap.add_argument("--beta", type=float, default=4.0); ap.add_argument("--dt", type=float, default=0.125)
    ap.add_argument("--tam", type=float, default=0.0); ap.add_argument("--t1", type=float, default=0.0)
    ap.add_argument("--tp", type=float, default=0.0, help="isotropic (spin-independent) NNN t'")
    ap.add_argument("--nwarm", type=int, default=300); ap.add_argument("--nmeas", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1); ap.add_argument("--ed", action="store_true")
    ap.add_argument("--chi", action="store_true", help="also measure tau-integrated pair susceptibility")
    ap.add_argument("--kres", action="store_true", help="k-resolve: emit maxk/k0/r0/rgt (unified_scan reduction)")
    a = ap.parse_args()
    q = DQMC(a.lx, a.ly, a.U, a.mu, a.beta, a.dt, a.tam, a.t1, a.seed, tp=a.tp)
    r = q.run(a.nwarm, a.nmeas, chi=a.chi, kres=a.kres)
    print(f"# DQMC {a.lx}x{a.ly} U={a.U} mu={a.mu} beta={a.beta} (NT={q.NT}) tam={a.tam} t1={a.t1} tp={a.tp}")
    print(f"  <sign> = {r['sign']:.4f}")
    print(f"  density = {r['dens']:.5f}")
    print(f"  energy(K+U) = {r['energy']:.5f}")
    print(f"  S_d (q=0) = {r['Sd']:.4f}   S_s = {r['Ss']:.4f}")
    if a.chi:
        print(f"  chi_d (q=0) = {r['chid']:.4f}   chi_s = {r['chis']:.4f}")
    if a.kres:
        cd = r["corr_d"]; print(f"  corr_d (full) maxk={cd['maxk']:.4f} k0={cd['k0']:.4f} r0={cd['r0']:.4f} rgt={cd['rgt']:.4f}")
        cv = r["corrV_d"]; print(f"  corrV_d (vtx) maxk={cv['maxk']:.4f} k0={cv['k0']:.4f} r0={cv['r0']:.4f} rgt={cv['rgt']:.4f}")
        if a.chi:
            xd = r["susc_d"]; print(f"  susc_d (full) maxk={xd['maxk']:.4f} k0={xd['k0']:.4f} r0={xd['r0']:.4f} rgt={xd['rgt']:.4f}")
            xv = r["suscV_d"]; print(f"  suscV_d (vtx) maxk={xv['maxk']:.4f} k0={xv['k0']:.4f} r0={xv['r0']:.4f} rgt={xv['rgt']:.4f}")
    if a.ed:
        e, d, Sd, Ss, chid, chis = ed_finite_T(a.lx, a.ly, a.U, a.mu, a.beta, a.tam, a.t1)
        print(f"  ED  density = {d:.5f}   energy = {e:.5f}   S_d = {Sd:.4f}   S_s = {Ss:.4f}")
        print(f"  ED  chi_d = {chid:.4f}   chi_s = {chis:.4f}")


if __name__ == "__main__":
    main()
