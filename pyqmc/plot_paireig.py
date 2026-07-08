import numpy as np, os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
root="/Users/macbook/Documents/Project/codex_CPQMC"
a=np.array([l.split(",") for l in open(os.path.join(root,"results/dqmc_scan/paireig_8x8_doped_b2.csv")) if l[0].isdigit()],float)
tam=a[:,0]; lam=a[:,1]; od=a[:,2]; os_=a[:,3]
fig,ax=plt.subplots(1,2,figsize=(12,5))
ax[0].plot(tam,lam,"C3-o",lw=2,ms=8,label=r"$\lambda_{\max}$ (leading pair eigenvalue)")
ax[0].axhline(1.0,ls="--",color="k",alpha=.4); ax[0].text(0.2,0.5,"instability at $\\lambda=1$",fontsize=9)
ax[0].set_ylim(0,1.05); ax[0].set_ylabel(r"$\lambda_{\max}$"); ax[0].set_title(r"leading pair eigenvalue $\approx0.002$ (no instability), flat in $t_{am}$")
ax[1].plot(tam,os_,"C0-o",lw=2,ms=8,label="extended-s")
ax[1].plot(tam,od,"C1-s",lw=2,ms=8,label="d-wave")
ax[1].set_ylim(0,1); ax[1].set_ylabel("leading eigenvector symmetry overlap")
ax[1].set_title("leading channel = extended-s, NOT d-wave")
for x in ax: x.set_xlabel(r"$t_{am}$ (altermagnet)",fontsize=12); x.grid(alpha=.3); x.legend(fontsize=11,frameon=False)
fig.suptitle(r"8$\times$8 pairing-eigenvalue (doped $n{\approx}0.88$, $\beta{=}2$, $\langle$sign$\rangle{=}1$): the direct test finds NO leading d-wave pairing, no $t_{am}$ enhancement",fontsize=12)
fig.tight_layout(rect=[0,0,1,0.95]); out=os.path.join(root,"docs","paireig_8x8_doped.png")
fig.savefig(out,dpi=140); print("wrote",out)
