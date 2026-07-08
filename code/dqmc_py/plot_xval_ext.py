#!/usr/bin/env python3
"""Cross-validation extensions: (A) tam axis -- T=0 CPQMC vs finite-T DQMC, showing the
SUSCEPTIBILITY suppressed but the equal-time CORRELATION enhanced (opposite to t1); and
(B) DQMC beta-series at t1=0 vs t1=0.2, showing chi_d^maxk growing toward the T=0 limit,
the t1 enhancement persisting at low T, and the doping sign problem (t1 protects <sign>)."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))
D = os.path.join(root, "results/dqmc_scan")


def ld_dqmc(name):
    a = np.array([l.split(",") for l in open(os.path.join(D, name)) if l[0].isdigit()], float)
    return dict(k=a[:, 0], sc_maxk=a[:, 3], cr_maxk=a[:, 6])


def ld_cpqmc(name):
    a = np.array([l.split() for l in open(os.path.join(D, name)) if l[0].isdigit()], float)
    return dict(k=a[:, 0], sc_maxk=a[:, 2], cr_maxk=a[:, 5])


def ld_beta(name):
    a = np.array([l.split() for l in open(os.path.join(D, name)) if l[0].isdigit()], float)
    return a  # beta, t1, dens, sign, susc_maxk


# ---------- Figure A: tam cross-validation ----------
dq = ld_dqmc("dqmc_tam_4x4_n0875.csv"); cp = ld_cpqmc("cpqmc_tam_4x4_n0875.txt")
figA, ax = plt.subplots(1, 2, figsize=(11, 4.4))
ax[0].plot(cp["k"], cp["sc_maxk"] / cp["sc_maxk"][0], "C0-o", lw=2, ms=8, label="T=0 CPQMC")
ax[0].plot(dq["k"], dq["sc_maxk"] / dq["sc_maxk"][0], "C3--s", lw=2, ms=8, mfc="none", label=r"finite-$T$ DQMC ($\beta{=}4$)")
ax[0].axhline(1, color="0.7", lw=.8, ls=":"); ax[0].set_title(r"peak-$q$ SUSCEPTIBILITY $\chi_d^{\max k}(t_{am})$")
ax[0].set_ylabel(r"normalized $\chi_d^{\max k}$"); ax[0].annotate("both: $t_{am}$\nSUPPRESSES", (0.25, 0.7), fontsize=9, color="0.25")
ax[1].plot(cp["k"], cp["cr_maxk"] / cp["cr_maxk"][0], "C0-o", lw=2, ms=8, label="T=0 CPQMC")
ax[1].plot(dq["k"], dq["cr_maxk"] / dq["cr_maxk"][0], "C3--s", lw=2, ms=8, mfc="none", label=r"finite-$T$ DQMC")
ax[1].axhline(1, color="0.7", lw=.8, ls=":"); ax[1].set_title(r"peak-$q$ equal-time CORRELATION $C_d^{\max k}(t_{am})$")
ax[1].set_ylabel(r"normalized $C_d^{\max k}$"); ax[1].annotate("both: $t_{am}$\nENHANCES", (0.25, 1.12), fontsize=9, color="0.25")
for a in ax:
    a.set_xlabel(r"$t_{am}$  (NN altermagnet)"); a.grid(alpha=.3); a.legend(frameon=False)
figA.suptitle(r"Cross-validation, $t_{am}$ axis: susceptibility DOWN, correlation UP -- both methods agree "
              r"(4$\times$4, $U{=}4$, $n{\approx}0.88$)", fontsize=11)
figA.tight_layout(rect=[0, 0, 1, 0.95])
outA = os.path.join(root, "docs", "xval_dqmc_cpqmc_tam.png"); figA.savefig(outA, dpi=130); print("wrote", outA)

# ---------- Figure B: beta-series ----------
b = ld_beta("dqmc_betaseries_4x4.txt")
b0 = b[b[:, 1] == 0.0]; b2 = b[b[:, 1] == 0.2]
figB, ax = plt.subplots(1, 2, figsize=(11, 4.4))
# left: susc vs beta (approach to T=0), marker size ~ sign reliability
ax[0].plot(b0[:, 0], b0[:, 4], "C7-o", lw=2, ms=7, label=r"$t_1=0$")
ax[0].plot(b2[:, 0], b2[:, 4], "C2-s", lw=2, ms=7, label=r"$t_1=0.2$")
ax[0].set_xlabel(r"$\beta = 1/T$"); ax[0].set_ylabel(r"$\chi_d^{\max k}$ (DQMC)")
ax[0].set_title(r"peak-$q$ $\chi_d$ grows toward $T{=}0$; $t_1$ stays above")
ax[0].annotate(r"$T\to 0$", (5.6, b0[:, 4].max() * 0.6), fontsize=10, color="0.4")
ax[0].grid(alpha=.3); ax[0].legend(frameon=False)
# right: sign vs beta (the doping sign problem; t1 helps)
ax[1].plot(b0[:, 0], b0[:, 3], "C7-o", lw=2, ms=7, label=r"$t_1=0$")
ax[1].plot(b2[:, 0], b2[:, 3], "C2-s", lw=2, ms=7, label=r"$t_1=0.2$")
ax[1].axhline(0.5, color="C3", lw=.8, ls=":"); ax[1].set_ylim(0, 1.03)
ax[1].set_xlabel(r"$\beta = 1/T$"); ax[1].set_ylabel(r"$\langle$sign$\rangle$")
ax[1].set_title(r"doping sign problem: $t_1$ protects $\langle$sign$\rangle$")
ax[1].annotate("shaded $\\beta$\nsign-limited", (4.8, 0.2), fontsize=8.5, color="C3")
ax[1].grid(alpha=.3); ax[1].legend(frameon=False)
figB.suptitle(r"finite-$T$ DQMC $\beta$-series (4$\times$4, $U{=}4$, $n{\approx}0.89$): approach to $T{=}0$ "
              r"and the doping sign problem", fontsize=11)
figB.tight_layout(rect=[0, 0, 1, 0.95])
outB = os.path.join(root, "docs", "dqmc_betaseries_4x4.png"); figB.savefig(outB, dpi=130); print("wrote", outB)
