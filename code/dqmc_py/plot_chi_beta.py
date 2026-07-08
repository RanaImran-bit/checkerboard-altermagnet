#!/usr/bin/env python3
"""Phase-1 capstone figure: DQMC tau-integrated d-wave pair susceptibility chi_d(beta)
vs finite-T ED (Kubo) on 2x2 U=4 half-filling, dt=0.03125. Shows the DQMC tracking ED
across temperature and growing toward the T=0 ground state as beta -> inf -- the
finite-T-DQMC -> T=0-CPQMC bridge. Numbers from the gated run (STATUS.md)."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

beta = np.array([1, 2, 3, 4, 5])
dqmc_chid = np.array([6.7096, 9.6259, 10.9319, 13.5195, 14.9997])
ed_chid   = np.array([6.5529, 10.1653, 11.8728, 13.2691, 14.7638])
dqmc_chis = np.array([3.6366, 3.7929, 3.8254, 3.8518, 3.85])   # beta=5 chis re-stated (parse glitch in scan)
ed_chis   = np.array([3.5809, 3.7643, 3.7672, 3.7624, 3.7579])

fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
ax[0].plot(beta, ed_chid, "k-o", lw=2, ms=6, label="finite-T ED (Kubo)")
ax[0].plot(beta, dqmc_chid, "C3--s", lw=2, ms=7, mfc="none", label="DQMC (this work)")
ax[0].set_title(r"$d$-wave pair susceptibility $\chi_d(\beta)$")
ax[1].plot(beta, ed_chis, "k-o", lw=2, ms=6, label="finite-T ED (Kubo)")
ax[1].plot(beta, dqmc_chis, "C0--s", lw=2, ms=7, mfc="none", label="DQMC (this work)")
ax[1].set_title(r"$s$-wave pair susceptibility $\chi_s(\beta)$")
for a in ax:
    a.set_xlabel(r"$\beta = 1/T$"); a.set_ylabel(r"$\chi$ (q=0)")
    a.legend(frameon=False); a.grid(alpha=0.3)
fig.suptitle(r"DQMC $\to$ ground state as $\beta\to\infty$  (2$\times$2 Hubbard, $U=4$, half-filling, $dt=0.03125$)",
             fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.96])
out = os.path.join(os.path.dirname(__file__), "..", "..", "docs", "dqmc_chi_beta.png")
fig.savefig(out, dpi=130)
print("wrote", os.path.normpath(out))
