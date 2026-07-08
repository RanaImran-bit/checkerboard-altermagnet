#!/usr/bin/env python3
"""Figure: finite-T DQMC maxk d-wave pair susceptibility vs the two altermagnet knobs
(tam, t1) at fixed doped filling (4x4, U=4, beta=4, n~0.94, sign>0.9). Independent-method
support for the manuscript: anisotropic t' (t1) ENHANCES the dynamic d-wave susceptibility
peak, while tam suppresses it and shifts the peak off q=0."""
import os, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(__file__)
root = os.path.normpath(os.path.join(here, "..", ".."))


def load(name):
    rows = []
    with open(os.path.join(root, "results", "dqmc_scan", name)) as f:
        for line in f:
            if line.startswith("#") or line.startswith("knob") or not line.strip():
                continue
            rows.append([float(x) for x in line.split(",")])
    a = np.array(rows)
    # cols: knob,dens,sign,sc_maxk,sc_k0,sc_r0,cr_maxk,cr_k0,cr_r0
    return dict(k=a[:, 0], dens=a[:, 1], sign=a[:, 2], sc_maxk=a[:, 3], sc_k0=a[:, 4],
               sc_r0=a[:, 5], cr_maxk=a[:, 6])


tam = load("dqmc_tam_4x4_b4.csv"); t1 = load("dqmc_t1_4x4_b4.csv")
fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.3), sharey=False)

for a, d, knob, col in ((ax[0], t1, r"$t_1$  (anisotropic NNN $t'$)", "C2"),
                        (ax[1], tam, r"$t_{am}$  (NN altermagnet)", "C3")):
    a.plot(d["k"], d["sc_maxk"], col + "-o", lw=2, ms=7, label=r"susc. peak-$q$  $\chi_d^{\max k}$")
    a.plot(d["k"], d["sc_k0"], col + "--s", lw=1.5, ms=6, mfc="none", label=r"susc. $q{=}0$  $\chi_d(0)$")
    a.plot(d["k"], d["cr_maxk"], "k:^", lw=1.3, ms=6, mfc="none", label=r"equal-time corr. peak-$q$")
    a.set_xlabel(knob); a.set_ylabel(r"$d$-wave pair amplitude")
    a.grid(alpha=0.3); a.legend(frameon=False, fontsize=8.5, loc="best")
    a2 = a.twinx(); a2.plot(d["k"], d["sign"], color="0.6", lw=1, marker=".")
    a2.set_ylabel(r"$\langle$sign$\rangle$", color="0.5"); a2.set_ylim(0.5, 1.02)
    a2.tick_params(axis="y", colors="0.5")

ax[0].set_title(r"$t_1$ ENHANCES peak-$q$ $\chi_d$ (susceptibility)")
ax[1].set_title(r"$t_{am}$ suppresses $\chi_d(0)$, peak shifts off $q{=}0$")
fig.suptitle(r"Finite-$T$ DQMC: $d$-wave pair susceptibility vs altermagnet knobs "
             r"(4$\times$4, $U{=}4$, $\beta{=}4$, $n{\approx}0.94$)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.95])
out = os.path.join(root, "docs", "dqmc_knob_scan.png")
fig.savefig(out, dpi=130)
print("wrote", out)
