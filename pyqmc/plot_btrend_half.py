import numpy as np, os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
root="/Users/macbook/Documents/Project/codex_CPQMC"
a=np.array([l.split(",") for l in open(os.path.join(root,"results/dqmc_scan/dqmc_btrend_half_8x8.csv")) if l[0].isdigit()],float)
# tam,t1,beta,corrF_maxk,corrF_k0,corrV_maxk,corrV_k0,dens,sign
pts=[(0.0,0.0,"baseline","C0"),(0.2,0.0,r"$t_{am}=0.2$","C1"),(0.4,0.0,r"$t_{am}=0.4$","C3"),(0.0,0.4,r"$t'=0.4$","C2")]
fig,ax=plt.subplots(1,2,figsize=(13,5.2))
for tam,t1,lab,c in pts:
    m=a[(a[:,0]==tam)&(a[:,1]==t1)]; b=m[:,2]; cV=m[:,5]; sg=m[:,8]
    ax[0].plot(b,cV,"-o",color=c,lw=2,ms=7,label=lab)
    # mark unreliable-sign points (sign<0.6) with open/red ring
    bad=sg<0.6
    if bad.any(): ax[0].plot(b[bad],cV[bad],"x",color="k",ms=11,mew=2)
    ax[1].plot(b,sg,"-o",color=c,lw=2,ms=7,label=lab)
ax[0].set_title(r"$d$-wave VERTEX corr $C_d^{V,\max q}$ vs $\beta$ (x = sign<0.6, untrustworthy)",fontsize=12)
ax[1].set_title(r"$\langle$sign$\rangle$ vs $\beta$ -- $t_{am}$ BREAKS the half-filling sign-free property",fontsize=12)
ax[1].axhline(0.6,ls="--",color="r",alpha=.5,lw=1)
for x in ax: x.set_xlabel(r"$\beta$",fontsize=13); x.grid(alpha=.3); x.legend(fontsize=11,frameon=False)
fig.suptitle("DQMC 8x8 HALF-FILLING (n=1) beta-trend: at T->0, tam SUPPRESSES d-wave vertex (and kills the sign)",fontsize=13)
fig.tight_layout(rect=[0,0,1,0.95]); out=os.path.join(root,"docs","dqmc_btrend_half_8x8.png")
fig.savefig(out,dpi=140); print("wrote",out)
