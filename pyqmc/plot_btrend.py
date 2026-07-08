#!/usr/bin/env python3
"""DQMC 8x8 d-wave pairing vs beta (equal-time), with the CPQMC T=0 reference, to check
the beta->inf convergence. Panels: corrF (full, should -> CPQMC T=0), corrV (vertex),
and <sign> vs beta (to flag where large-beta data is untrustworthy)."""
import os, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
root="/Users/macbook/Documents/Project/codex_CPQMC"; D=os.path.join(root,"results/dqmc_scan")
a=np.array([l.split(",") for l in open(os.path.join(D,"dqmc_btrend_8x8.csv")) if l[0].isdigit()],float)
# cols: tam,t1,beta,corrF_maxk,corrF_k0,corrV_maxk,corrV_k0,dens,sign
cp=np.array([l.split(",") for l in open(os.path.join(D,"cpqmc_vertex2d_8x8.csv")) if l[0].isdigit()],float)
def cpval(tam,t1,col):  # CPQMC T=0 reference (col 16=corrF_maxk,12=corrV_maxk 0-idx)
    m=cp[(cp[:,0]==tam)&(cp[:,1]==t1)]; return m[0,col] if len(m) else np.nan
pts=[(0.0,0.0,"baseline"),(0.2,0.0,r"$t_{am}{=}0.2$"),(0.4,0.0,r"$t_{am}{=}0.4$"),(0.0,0.4,r"$t'{=}0.4$")]
fig,ax=plt.subplots(1,3,figsize=(16,5))
for tam,t1,lab in pts:
    m=a[(a[:,0]==tam)&(a[:,1]==t1)]; b=m[:,2]
    ax[0].plot(b,m[:,3],"-o",lw=2,ms=7,label=lab)
    ax[1].plot(b,m[:,5],"-o",lw=2,ms=7,label=lab)
    ax[2].plot(b,m[:,8],"-o",lw=2,ms=7,label=lab)
    ax[0].axhline(cpval(tam,t1,16),ls=":",color=ax[0].lines[-1].get_color(),lw=1.2)
    ax[1].axhline(cpval(tam,t1,12),ls=":",color=ax[1].lines[-1].get_color(),lw=1.2)
ax[0].set_title(r"FULL corr $C_d^{\max q}$ (dotted = CPQMC T=0)",fontsize=13)
ax[1].set_title(r"VERTEX corr $C_d^{V,\max q}$ (dotted = CPQMC T=0)",fontsize=13)
ax[2].set_title(r"$\langle$sign$\rangle$ vs $\beta$ (trust where ~1)",fontsize=13)
for x in ax: x.set_xlabel(r"$\beta$",fontsize=13); x.grid(alpha=.3); x.legend(fontsize=10,frameon=False)
ax[2].axhline(0.5,ls="--",color="r",lw=1,alpha=.5)
fig.suptitle("DQMC 8x8 (n~0.97) d-wave pairing vs beta -> does it converge to CPQMC T=0?",fontsize=14)
fig.tight_layout(rect=[0,0,1,0.95]); out=os.path.join(root,"docs","dqmc_btrend_8x8.png")
fig.savefig(out,dpi=140); print("wrote",out)
