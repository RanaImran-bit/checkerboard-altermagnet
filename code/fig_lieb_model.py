"""The Lieb lattice as a second emergent-altermagnet model.

Three sites per cell: A on the square-lattice corners, B on the horizontal bond
midpoints, C on the vertical ones. Nearest-neighbour hopping connects A-B and
A-C only, and is SPIN-INDEPENDENT, so at U = 0 the two spin species are exactly
degenerate -- there is no splitting to inherit. B and C are related by a 90
degree ROTATION (x <-> y), which is the altermagnet criterion.

    H0(k) = [[0,   a,   b ],      a = 2t cos(kx/2)
             [a,   0,   0 ],      b = 2t cos(ky/2)
             [b,   0,   0 ]]

Altermagnetic order puts OPPOSITE moments on B and C (A stays unpolarised):

    H_sigma(k) = H0(k) + sigma * m * diag(0, +1, -1)

Swapping kx <-> ky swaps a <-> b, which exchanges the roles of B and C and is
therefore equivalent to flipping sigma. So the splitting is ODD under the
diagonal mirror -- d_x2-y2 symmetry. That is the OPPOSITE d channel from our
checkerboard, which came out d_xy. Two emergent-AM models, the two different
d-wave symmetries, same mechanism.
"""
import numpy as np, matplotlib as mpl, matplotlib.pyplot as plt

mpl.rcParams.update({
    'font.family':'serif','font.serif':['DejaVu Serif'],
    'mathtext.fontset':'dejavuserif','axes.linewidth':1.2,
    'axes.labelsize':15,'xtick.labelsize':12,'ytick.labelsize':12})

t=-1.0
def Hk(kx,ky,m,sigma):
    a=2*t*np.cos(kx/2); b=2*t*np.cos(ky/2)
    H=np.array([[0,a,b],[a,sigma*m,0],[b,0,-sigma*m]],dtype=float)
    return np.linalg.eigvalsh(H)

fig=plt.figure(figsize=(17.5,9.6),facecolor='white')
gs=fig.add_gridspec(2,3,hspace=0.30,wspace=0.28)

# ---------------- (a) lattice schematic ----------------
a0=fig.add_subplot(gs[0,0])
for X in range(3):
    for Y in range(3):
        a0.plot([X,X+1],[Y,Y],'k-',lw=1.0,zorder=1)
        a0.plot([X,X],[Y,Y+1],'k-',lw=1.0,zorder=1)
for X in range(4):
    for Y in range(4):
        a0.plot(X,Y,'o',ms=13,mfc='0.75',mec='k',mew=1.2,zorder=3)          # A
for X in range(3):
    for Y in range(4):
        a0.plot(X+0.5,Y,'o',ms=13,mfc='tab:red',mec='k',mew=1.2,zorder=3)   # B
for X in range(4):
    for Y in range(3):
        a0.plot(X,Y+0.5,'o',ms=13,mfc='tab:blue',mec='k',mew=1.2,zorder=3)  # C
a0.plot([],[],'o',ms=11,mfc='0.75',mec='k',label='A  (corner, unpolarised)')
a0.plot([],[],'o',ms=11,mfc='tab:red',mec='k',label=r'B  ($x$-bond, $+m$)')
a0.plot([],[],'o',ms=11,mfc='tab:blue',mec='k',label=r'C  ($y$-bond, $-m$)')
a0.legend(fontsize=10,loc='upper center',bbox_to_anchor=(0.5,-0.06),frameon=False)
a0.set_xlim(-0.6,3.6); a0.set_ylim(-0.6,3.6); a0.set_aspect('equal'); a0.axis('off')
a0.set_title('(a)  Lieb lattice: B and C related by $90^\\circ$ rotation',fontsize=14)

# ---------------- (b,c) bands along Gamma-X-M-Gamma ----------------
def path(n=180):
    G=(0,0); X=(np.pi,0); M=(np.pi,np.pi)
    segs=[(G,X),(X,M),(M,G)]; ks=[]; ticks=[0]
    for s,e in segs:
        for i in range(n):
            ks.append((s[0]+(e[0]-s[0])*i/n, s[1]+(e[1]-s[1])*i/n))
        ticks.append(len(ks))
    return ks,ticks
ks,ticks=path()
for col,m,lab in [(1,0.0,'(b)  $m=0$  (i.e. $U=0$): spin degenerate'),
                  (2,0.8,'(c)  $m=0.8$: spin split')]:
    ax=fig.add_subplot(gs[0,col])
    for sig,c,ls,nm in [(+1,'tab:red','-',r'$\uparrow$'),(-1,'tab:blue','--',r'$\downarrow$')]:
        E=np.array([Hk(kx,ky,m,sig) for kx,ky in ks])
        for band in range(3):
            ax.plot(E[:,band],color=c,ls=ls,lw=1.8,label=nm if band==0 else None)
    ax.set_xticks(ticks); ax.set_xticklabels([r'$\Gamma$','X','M',r'$\Gamma$'])
    ax.set_ylabel(r'$E/t$'); ax.set_title(lab,fontsize=14)
    ax.axhline(0,color='0.6',lw=0.8,ls=':'); ax.legend(fontsize=12,frameon=False)
    ax.grid(alpha=0.15)

# ---------------- (d,e) splitting over the BZ ----------------
N=240
kx=np.linspace(-np.pi,np.pi,N); ky=np.linspace(-np.pi,np.pi,N)
KX,KY=np.meshgrid(kx,ky)
for col,m in [(0,0.0),(1,0.8)]:
    ax=fig.add_subplot(gs[1,col])
    S=np.zeros_like(KX)
    for i in range(N):
        for j in range(N):
            up=Hk(KX[i,j],KY[i,j],m,+1); dn=Hk(KX[i,j],KY[i,j],m,-1)
            S[i,j]=up[1]-dn[1]                      # middle band
    v=max(np.abs(S).max(),1e-12)
    im=ax.contourf(KX,KY,S,levels=60,cmap='RdBu_r',vmin=-v,vmax=v)
    fig.colorbar(im,ax=ax,pad=0.02)
    ax.set_aspect('equal'); ax.set_xlabel(r'$k_x$'); ax.set_ylabel(r'$k_y$')
    ax.set_xticks([-np.pi,0,np.pi]); ax.set_xticklabels([r'$-\pi$','0',r'$\pi$'])
    ax.set_yticks([-np.pi,0,np.pi]); ax.set_yticklabels([r'$-\pi$','0',r'$\pi$'])
    ax.set_title(rf'({"de"[col]})  $E_\uparrow-E_\downarrow$,  $m={m:g}$   '
                 rf'(max $={np.abs(S).max():.3f}$)',fontsize=13)
    print(f"m={m}: max |E_up - E_dn| = {np.abs(S).max():.6f}")

# ---------------- (f) symmetry test of the splitting ----------------
ax=fig.add_subplot(gs[1,2])
m=0.8
S=np.zeros_like(KX)
for i in range(N):
    for j in range(N):
        S[i,j]=Hk(KX[i,j],KY[i,j],m,+1)[1]-Hk(KX[i,j],KY[i,j],m,-1)[1]
rot=np.rot90(S)                      # 90 degree rotation
mir=S.T                              # mirror about kx = ky
A_odd=np.abs(S+rot).sum()/np.abs(S).sum()
M_odd=np.abs(S+mir).sum()/np.abs(S).sum()
ax.axis('off')
ax.text(0.02,0.86,'(f)  symmetry of the splitting',fontsize=14,transform=ax.transAxes)
ax.text(0.05,0.66,rf'$A_{{\rm odd}} = {A_odd:.3f}$   (0 = C4-odd $\Rightarrow$ $d$-wave)',
        fontsize=13,transform=ax.transAxes)
ax.text(0.05,0.50,rf'$M_{{\rm odd}} = {M_odd:.3f}$   (0 = $d_{{x^2-y^2}}$,  2 = $d_{{xy}}$)',
        fontsize=13,transform=ax.transAxes)
ax.text(0.05,0.28,'Lieb $\\Rightarrow$ $d_{x^2-y^2}$\ncheckerboard $\\Rightarrow$ $d_{xy}$\n'
        'same mechanism, opposite $d$ channel',fontsize=13,transform=ax.transAxes,color='tab:red')
print(f"A_odd={A_odd:.4f}  M_odd={M_odd:.4f}")
plt.savefig('/Users/liujiaxin/Desktop/checkerboard-altermagnet/figs/fig_lieb_model.png',
            dpi=300,bbox_inches='tight',facecolor='white')
