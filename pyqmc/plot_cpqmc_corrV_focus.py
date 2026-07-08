#!/usr/bin/env python3
"""Clear appeal heatmap: T=0 CPQMC equal-time d-wave pairing VERTEX over (t', tam).
The USABLE enhancement is tam (altermagnet), NOT t': tam raises the peak-q vertex and
drives q=0 negative (pairing -> finite q); t' suppresses everywhere."""
import os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 13})
root="/Users/macbook/Documents/Project/codex_CPQMC"
csv=sys.argv[1]; tag=csv.replace("cpqmc_vertex2d_","").replace(".csv","")
a=np.array([l.split(",") for l in open(os.path.join(root,"results/dqmc_scan",csv)) if l[0].isdigit()],float)
tam=np.unique(a[:,0]); t1=np.unique(a[:,1]); nt,n1=len(tam),len(t1); dens=a[0,2]
ext=[t1.min()-.05,t1.max()+.05,tam.min()-.05,tam.max()+.05]
fig,ax=plt.subplots(1,2,figsize=(14.5,6.2))
for k,(col,ttl,cm) in enumerate([(12,r"VERTEX pairing at peak $q$   $C_d^{\,\mathrm{max}\,q}$",  "magma"),
                                 (13,r"VERTEX pairing at $q=0$   $C_d^{\,q=0}$","magma")]):
    M=a[:,col].reshape(nt,n1)
    im=ax[k].imshow(M,origin="lower",extent=ext,cmap=cm,aspect="auto",interpolation="bilinear")
    cs=ax[k].contour(t1,tam,M,5,colors="w",linewidths=0.6,alpha=0.5)
    for i,tv in enumerate(tam):
        for j,t1v in enumerate(t1):
            ax[k].text(t1v,tv,f"{M[i,j]:.0f}",ha="center",va="center",fontsize=12,fontweight="bold",
                       color="w" if M[i,j]<0.55*M.max() else "k")
    ax[k].set_xlabel(r"$t_1$   (anisotropic $t'$)  $\longrightarrow$ suppresses",fontsize=13)
    ax[k].set_ylabel(r"$t_{am}$   (altermagnet)  $\longrightarrow$",fontsize=13)
    ax[k].set_title(ttl,fontsize=14,pad=8)
    cb=fig.colorbar(im,ax=ax[k],fraction=.046,pad=.04)
    ax[k].annotate("",xy=(-0.02,0.42),xytext=(-0.02,0.0),
                   arrowprops=dict(arrowstyle="-|>",color="lime",lw=3),annotation_clip=False)
ax[0].text(-0.16,0.2,"tam ENHANCES",rotation=90,color="green",fontsize=13,fontweight="bold",
           va="center",ha="center",transform=ax[0].transData)
fig.suptitle(rf"T=0 CPQMC  {tag.split('_')[0]}  ($U{{=}}4$, $n{{=}}{dens:.3f}$):  the ALTERMAGNET ($t_{{am}}$) "
             rf"grows finite-$q$ $d$-wave pairing; $t'$ does not",fontsize=15)
fig.tight_layout(rect=[0,0,1,0.94])
out=os.path.join(root,"docs",f"cpqmc_corrV_focus_{tag}.png"); fig.savefig(out,dpi=150); print("wrote",out)
