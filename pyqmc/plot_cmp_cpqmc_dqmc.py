import numpy as np, os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
D="results/dqmc_scan"
def load(f): 
    return np.array([l.split(",") for l in open(os.path.join(D,f)) if l[0].isdigit()],float)
cp=load("cpqmc_vertex2d_8x8.csv"); d2=load("dqmc_vertex2d_8x8_b2.csv"); d4=load("dqmc_vertex2d_8x8_b4.csv")
def axis(a,fix_t1):  # return knob, corrV_maxk along tam (t1=0) or t1 (tam=0)
    if fix_t1: m=a[a[:,1]==0.0]; x=m[:,0]
    else: m=a[a[:,0]==0.0]; x=m[:,1]
    return x, m[:,12]   # corrV maxk
fig,ax=plt.subplots(1,2,figsize=(13,5))
for k,(fix,lab) in enumerate([(True,r"$t_{am}$ (altermagnet)"),(False,r"$t_1$ (anisotropic $t'$)")]):
    for a,nm,st in [(cp,"CPQMC T=0","C0-o"),(d2,r"DQMC $\beta{=}2$","C3--s"),(d4,r"DQMC $\beta{=}4$","C2:^")]:
        x,y=axis(a,fix); ax[k].plot(x,y,st,lw=2,ms=8,mfc="none" if st[0]!="C0" else None,label=nm)
    ax[k].set_xlabel(lab,fontsize=13); ax[k].set_ylabel(r"vertex $C_d$ at peak $q$",fontsize=13)
    ax[k].grid(alpha=.3); ax[k].legend(fontsize=11,frameon=False)
ax[0].set_title("along tam: all rise (enhance)",fontsize=13)
ax[1].set_title("along t': all fall (suppress)",fontsize=13)
fig.suptitle("CPQMC vs DQMC, 8x8 vertex pairing: TRENDS agree, magnitudes differ (T-dependent + noise)",fontsize=14)
fig.tight_layout(rect=[0,0,1,0.95])
out="docs/cmp_cpqmc_dqmc_8x8.png"; fig.savefig(out,dpi=140); print("wrote",out)
