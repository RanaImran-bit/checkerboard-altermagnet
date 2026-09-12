import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, matplotlib.pyplot as plt
D='/private/tmp/claude-501/-Users-liujiaxin-Desktop-Susceptibility-qmc-platform-master/c2fafe5f-4f29-467d-a412-3fdcefa75ac0/scratchpad'
d=pd.read_csv(f'{D}/pr/pr_U4_all.csv'); d['n']=2*d.nup/144
CH=[('son','on-site s','#7f7f7f'),('sext','extended s','#2ca02c'),
    ('d','$d_{x^2-y^2}$','#1f77b4'),('dxy','$d_{xy}$','#d62728')]
fig,ax=plt.subplots(1,2,figsize=(13.5,4.8),facecolor='white')
# left: local vs long-range magnitude, half filling, delta = 0
h=d[np.isclose(d.n,1.0)]
w=0.35; x=np.arange(4)
for k,dd in enumerate([0.0,0.4]):
    loc=[h[(h.channel==c)&(h.R==0)&np.isclose(h.delta,dd)].V.mean() for c,_,_ in CH]
    ax[0].bar(x+(k-0.5)*w,loc,w,label=rf'$\delta={dd:g}$',
              color=['#bbb','#C0392B'][k],edgecolor='k',linewidth=0.7)
ax[0].axhline(0,color='k',lw=1); ax[0].set_xticks(x)
ax[0].set_xticklabels([l for _,l,_ in CH],fontsize=13)
ax[0].set_ylabel(r'local  $V_\alpha(R=0)$',fontsize=15)
ax[0].set_title('(a)  local vertex',fontsize=16); ax[0].legend(fontsize=12)
for k,dd in enumerate([0.0,0.4]):
    lng=[h[(h.channel==c)&(h.R>2)&np.isclose(h.delta,dd)].V.mean() for c,_,_ in CH]
    ax[1].bar(x+(k-0.5)*w,lng,w,label=rf'$\delta={dd:g}$',
              color=['#bbb','#C0392B'][k],edgecolor='k',linewidth=0.7)
ax[1].axhline(0,color='k',lw=1); ax[1].set_xticks(x)
ax[1].set_xticklabels([l for _,l,_ in CH],fontsize=13)
ax[1].set_ylabel(r'long-range  $\overline{V_\alpha}(R>2)$',fontsize=15)
ax[1].set_title('(b)  long-range vertex',fontsize=16); ax[1].legend(fontsize=12)
for a in ax: a.tick_params(labelsize=12)
fig.suptitle(r'Anisotropy moves the long-range attraction from $d_{x^2-y^2}$ to $d_{xy}$'
             '\n'+r'$L=12$, $U=4$, half filling',fontsize=16)
plt.tight_layout(rect=[0,0,1,0.87])
plt.savefig(f'{D}/deck2/img/prbar.png',dpi=110,facecolor='white')
print('bar figure written')
