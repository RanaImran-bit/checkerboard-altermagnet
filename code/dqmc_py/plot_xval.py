#!/usr/bin/env python3
"""Phase-3 cross-validation figure: the maxk d-wave PAIR SUSCEPTIBILITY vs anisotropic
NNN t' (=t1), computed two independent ways on the SAME 4x4 Hubbard lattice at matched
filling n~0.88 -- T=0 constrained-path CPQMC and finite-T (beta=4) DQMC. Both show the
SAME enhancement with t1; the equal-time correlation moves oppositely. Absolute scales
differ (different tau-integration window + T), so the left panel normalizes each curve to
its t1=0 value to overlay the trend; the right panel shows the equal-time correlation."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))


def load_dqmc(name):
    rows = [l.split(",") for l in open(os.path.join(root, "results/dqmc_scan", name))
            if l[0].isdigit()]
    a = np.array(rows, float)
    return dict(k=a[:, 0], sc_maxk=a[:, 3], cr_maxk=a[:, 6])


def load_cpqmc(name):
    rows = [l.split() for l in open(os.path.join(root, "results/dqmc_scan", name))
            if l[0].isdigit()]
    a = np.array(rows, float)
    return dict(k=a[:, 0], sc_maxk=a[:, 2], cr_maxk=a[:, 5])


dq = load_dqmc("dqmc_t1_4x4_n0875.csv")
cp = load_cpqmc("cpqmc_t1_4x4_n0875.txt")

fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))

# left: susceptibility, normalized to t1=0 (trend overlay)
ax[0].plot(cp["k"], cp["sc_maxk"] / cp["sc_maxk"][0], "C0-o", lw=2, ms=8,
           label=r"$T{=}0$ CPQMC (constrained path)")
ax[0].plot(dq["k"], dq["sc_maxk"] / dq["sc_maxk"][0], "C3--s", lw=2, ms=8, mfc="none",
           label=r"finite-$T$ DQMC ($\beta{=}4$)")
ax[0].axhline(1.0, color="0.7", lw=0.8, ls=":")
ax[0].set_title(r"peak-$q$ $d$-wave SUSCEPTIBILITY  $\chi_d^{\max k}(t_1)/\chi_d^{\max k}(0)$")
ax[0].set_ylabel(r"normalized $\chi_d^{\max k}$")
ax[0].annotate("both methods:\nt' ENHANCES", (0.27, 1.05), fontsize=9, color="0.25")

# right: equal-time correlation, normalized
ax[1].plot(cp["k"], cp["cr_maxk"] / cp["cr_maxk"][0], "C0-o", lw=2, ms=8,
           label=r"$T{=}0$ CPQMC")
ax[1].plot(dq["k"], dq["cr_maxk"] / dq["cr_maxk"][0], "C3--s", lw=2, ms=8, mfc="none",
           label=r"finite-$T$ DQMC")
ax[1].axhline(1.0, color="0.7", lw=0.8, ls=":")
ax[1].set_title(r"peak-$q$ equal-time CORRELATION  $C_d^{\max k}(t_1)/C_d^{\max k}(0)$")
ax[1].set_ylabel(r"normalized $C_d^{\max k}$")
ax[1].annotate("both methods:\nt' decreases", (0.05, 0.86), fontsize=9, color="0.25")

for a in ax:
    a.set_xlabel(r"$t_1$  (anisotropic NNN $t'$)"); a.grid(alpha=0.3); a.legend(frameon=False)
fig.suptitle(r"Cross-validation: $d$-wave pairing vs $t'$ -- T=0 CPQMC and finite-$T$ DQMC agree "
             r"(4$\times$4, $U{=}4$, $n{\approx}0.88$)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.95])
out = os.path.join(root, "docs", "xval_dqmc_cpqmc_t1.png")
fig.savefig(out, dpi=130)
print("wrote", out)
