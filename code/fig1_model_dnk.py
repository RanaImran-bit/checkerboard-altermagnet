"""Fig 1: checkerboard model schematic (a) + mean-field Delta n(k) dxy maps (b,c,d).
Non-interacting/mean-field, no QMC data needed. Runs anywhere."""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.transforms import blended_transform_factory
from matplotlib.patches import FancyArrowPatch
from matplotlib.lines import Line2D
import os

plt.rcParams.update({
    'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans','sans-serif'],
    'mathtext.fontset':'dejavusans','axes.unicode_minus':False,'axes.linewidth':2.0,
    'xtick.major.width':2.0,'ytick.major.width':2.0,'xtick.major.size':7,'ytick.major.size':7,
})

def plot_model_combined(ax, fontsize=18):
    UP='#c0392b'; DN='#1a5276'; GREY='#9e9e9e'; TPLUS='#2a9d3a'; TMINUS='#e08e0b'
    n=4
    for i in range(n):
        for j in range(n):
            if i<n-1: ax.plot([i,i+1],[j,j],'-',color=GREY,lw=1.6,zorder=1)
            if j<n-1: ax.plot([i,i],[j,j+1],'-',color=GREY,lw=1.6,zorder=1)
    for i in range(n-1):
        for j in range(n-1):
            if (i+j)%2==0:
                ax.plot([i,i+1],[j,j+1],'-',color=TPLUS,lw=3.0,zorder=2)
                ax.plot([i+1,i],[j,j+1],'-',color=TPLUS,lw=3.0,zorder=2)
            else:
                ax.plot([i,i+1],[j,j+1],'--',color=TMINUS,lw=2.6,dashes=(5,2),zorder=2)
                ax.plot([i+1,i],[j,j+1],'--',color=TMINUS,lw=2.6,dashes=(5,2),zorder=2)
    def spin_arrow(x,y,up):
        c=UP if up else DN; dy=0.34 if up else -0.34; y0=y-dy/2
        ax.add_patch(FancyArrowPatch((x,y0),(x,y0+dy),arrowstyle='-|>',mutation_scale=15,lw=2.4,color=c,zorder=6))
    for i in range(n):
        for j in range(n):
            up=(i+j)%2==0; c=UP if up else DN
            ax.plot(i,j,'o',ms=22,mfc='white',mec=c,mew=2.2,zorder=5); spin_arrow(i,j,up)
    leg=[Line2D([0],[0],color=TPLUS,lw=3.0,label=r"$t'+\delta$"),
         Line2D([0],[0],color=TMINUS,lw=2.6,ls='--',label=r"$t'-\delta$"),
         Line2D([0],[0],color=GREY,lw=1.6,label=r"$t$")]
    ax.legend(handles=leg,loc='lower center',bbox_to_anchor=(0.5,-0.08),
              fontsize=fontsize-3,frameon=False,ncol=3,handlelength=1.8,columnspacing=1.2,handletextpad=0.5)
    ax.set_xlim(-0.7,n-0.3); ax.set_ylim(-0.9,n-0.3); ax.set_aspect('equal'); ax.axis('off')
    ax.text(0.0,1.0,'(a)',transform=ax.transAxes,fontsize=fontsize+6,fontweight='bold',va='top',ha='left')

def eps_up(kx,ky,delta,M,t=1.0,tp=-0.3):
    e0=-4*tp*np.cos(kx)*np.cos(ky); g=2*t*(np.cos(kx)+np.cos(ky)); m=4*delta*np.sin(kx)*np.sin(ky)-M
    return e0-np.sqrt(g**2+m**2)
def eps_dn(kx,ky,delta,M,t=1.0,tp=-0.3):
    e0=-4*tp*np.cos(kx)*np.cos(ky); g=2*t*(np.cos(kx)+np.cos(ky)); m=4*delta*np.sin(kx)*np.sin(ky)+M
    return e0-np.sqrt(g**2+m**2)
def fermi(E,mu,T): return 1.0/(1.0+np.exp((E-mu)/T))

N=500; k=np.linspace(-np.pi,np.pi,N); KX,KY=np.meshgrid(k,k)
M0=1.0; MU=-1.6; T=0.15                     # schematic moment, filling, broadening
cases=[(0.0,r'$\delta = 0$'),(0.2,r'$\delta = 0.2$'),(0.4,r'$\delta = 0.4$')]

fig=plt.figure(figsize=(23,7.4),facecolor='white')
gs=GridSpec(1,5,figure=fig,width_ratios=[1.25,1,1,1,0.06],wspace=0.16,left=0.05,right=0.95,top=0.86,bottom=0.20)
ax_lat=fig.add_subplot(gs[0,0]); axs=[fig.add_subplot(gs[0,c+1]) for c in range(3)]; cax=fig.add_subplot(gs[0,4])
plot_model_combined(ax_lat,fontsize=18)

TLBL=26; PL=22; fs_lbl=['(b)','(c)','(d)']; im=None
for ci,(d,ttl) in enumerate(cases):
    EU=eps_up(KX,KY,d,M0); ED=eps_dn(KX,KY,d,M0)
    dn=fermi(EU,MU,T)-fermi(ED,MU,T)
    ax=axs[ci]
    im=ax.imshow(dn,extent=[-np.pi,np.pi,-np.pi,np.pi],origin='lower',cmap='RdBu_r',vmin=-1,vmax=1,aspect='equal',zorder=1)
    for px,py,nm,dx,dy in [(0,0,r'$\Gamma$',0.12,0.10),(np.pi,np.pi,r'$M$',-0.16,-0.14),(np.pi,0,r'$X$',-0.12,0.10)]:
        ax.scatter(px,py,s=30,c='black',zorder=5)
        ax.text(px+dx*np.pi,py+dy*np.pi,nm,fontsize=20,fontweight='bold',ha='center',va='center',zorder=6)
    for sp in ax.spines.values(): sp.set_color('black'); sp.set_linewidth(2.0)
    ax.set_xticks([-np.pi,0,np.pi]); ax.set_yticks([-np.pi,0,np.pi]); ax.set_xticklabels([]); ax.set_yticklabels([])
    tx=blended_transform_factory(ax.transData,ax.transAxes)
    for xv,lab in [(-np.pi,r'$-\pi$'),(0,r'$0$'),(np.pi,r'$\pi$')]:
        ax.text(xv,-0.10,lab,transform=tx,ha='center',va='top',fontsize=TLBL,fontweight='bold',clip_on=False)
    if ci==0:
        ty=blended_transform_factory(ax.transAxes,ax.transData)
        for yv,lab,va_a in [(-np.pi,r'$-\pi$','bottom'),(0,r'$0$','center'),(np.pi,r'$\pi$','top')]:
            ax.text(-0.08,yv,lab,transform=ty,ha='right',va=va_a,fontsize=TLBL,fontweight='bold',clip_on=False)
        ax.text(-0.24,0.5,r'$k_y$',transform=ax.transAxes,ha='center',va='center',fontsize=TLBL+4,fontweight='bold',rotation=90)
    if ci==1:
        ax.text(0.5,-0.19,r'$k_x$',transform=ax.transAxes,ha='center',va='top',fontsize=TLBL+4,fontweight='bold')
    ax.set_title(ttl,fontsize=24,fontweight='bold',pad=10)
    ax.text(0.04,0.96,fs_lbl[ci],transform=ax.transAxes,fontsize=PL,fontweight='bold',va='top',color='black',
            bbox=dict(boxstyle='round,pad=0.15',fc='white',ec='none',alpha=0.7))
cb=fig.colorbar(im,cax=cax); cb.set_ticks([-1,0,1])
cb.set_label(r'$\Delta n(\mathbf{k}) = n_\uparrow(\mathbf{k})-n_\downarrow(\mathbf{k})$',fontsize=17)
cax.text(0.5,1.03,r'spin $\uparrow$',transform=cax.transAxes,ha='center',va='bottom',fontsize=13,color='#7a1620',fontweight='bold')
cax.text(0.5,-0.03,r'spin $\downarrow$',transform=cax.transAxes,ha='center',va='top',fontsize=13,color='#123a55',fontweight='bold')

OUT=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','figures','fig1_model_dnk.png')
fig.savefig(OUT,dpi=300,bbox_inches='tight',facecolor='white'); print("saved",OUT)
