"""Panel (b) at L = 8, 10 and 12: is the reversal at delta ~ 0.3 finite size?

x = Delta_tot/N, y = anisotropy delta, colour = U. One row per delta; the colour
order along that row gives the SIGN of dDelta_tot/dU:

    dark on the left, bright on the right  ->  Delta_tot GROWS with U
    bright on the left, dark on the right  ->  Delta_tot FALLS with U

The reversal claim is that the sign flips between delta < 0.3 and delta > 0.3.
If that flip survives at all three L it is physical; if it moves or vanishes it
was a shell artifact.

Per site (Eq. 3 sum / L^2), the only form comparable across L. All three are
Fortran runs at the production settings, so nothing here is affected by the
beta = 3 problem in the Python chi drivers.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt

mpl.rcParams.update({
    "font.family":"serif","font.serif":["DejaVu Serif"],
    "mathtext.fontset":"dejavuserif","axes.linewidth":1.1,
    "xtick.direction":"in","ytick.direction":"in",
    "xtick.top":True,"ytick.right":True,
    "axes.labelsize":16,"xtick.labelsize":13,"ytick.labelsize":13})

D="/Users/liujiaxin/Desktop/checkerboard-altermagnet/data"
d=pd.read_csv(f"{D}/fortran_L8_L10_L12_dedup.csv")
d=d[np.isclose(d.n,1.0)&(d.U>=2)]
LS=sorted(d.L.unique())
norm=plt.Normalize(d.U.min(),d.U.max())

fig,ax=plt.subplots(1,len(LS),figsize=(6.0*len(LS),5.2),sharey=True,facecolor="white")
for i,L in enumerate(LS):
    a=ax[i]; s=d[d.L==L]
    sc=a.scatter(s.dtot_N,s.delta,c=s.U,cmap="viridis",norm=norm,s=300,
                 edgecolor="0.2",linewidth=0.7,zorder=3)
    a.set_xlabel(r"$\Delta_{\mathrm{tot}}/N$")
    if i==0: a.set_ylabel(r"anisotropy  $\delta$")
    a.set_title(rf"$L={L}$",fontsize=17)
    a.grid(alpha=0.18,lw=0.7); a.set_axisbelow(True)
cb=fig.colorbar(sc,ax=ax,pad=0.015,fraction=0.018)
cb.set_label(r"$U/t$",fontsize=16)
plt.savefig(f"{D}/../figs/fig_panelb_3L.png",dpi=600,bbox_inches="tight",facecolor="white")

print("sign of dDelta_tot/dU at each delta  (corr of Delta_tot with U, fixed delta)\n")
print(f"{'delta':>7}" + "".join(f"{'L='+str(L):>12s}" for L in LS) + "   verdict")
for dl in sorted(d.delta.unique()):
    row=f"{dl:7.1f}"; rs=[]
    for L in LS:
        g=d[(d.L==L)&(d.delta==dl)]
        r=np.corrcoef(g.U,g.dtot_N)[0,1]; rs.append(r)
        row+=f"{r:+12.2f}"
    sg={np.sign(round(r,2)) for r in rs}
    row+="   " + ("consistent" if len(sg)==1 else "DIFFERS with L")
    print(row)
