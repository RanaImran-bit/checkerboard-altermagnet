import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, matplotlib.pyplot as plt
D='/private/tmp/claude-501/-Users-liujiaxin-Desktop-Susceptibility-qmc-platform-master/c2fafe5f-4f29-467d-a412-3fdcefa75ac0/scratchpad/pr'
CH=[('son',r'on-site $s$','#7f7f7f','o'),('sext',r'extended $s$','#2ca02c','s'),
    ('d',r'$d_{x^2-y^2}$','#1f77b4','^'),('dxy',r'$d_{xy}$','#d62728','D')]
d=pd.read_csv(f'{D}/pr_U4_all.csv'); d['n']=2*d.nup/144
h=d[np.isclose(d.n,1.0)]

fig,ax=plt.subplots(1,3,figsize=(18,5.0),facecolor='white')

# (a) local, R = 0
for c,lab,col,mk in CH:
    s=h[(h.channel==c)&(h.R==0)].groupby('delta').V.agg(['mean','sem'])
    ax[0].errorbar(s.index,s['mean'],s['sem'],marker=mk,ms=8,lw=2,color=col,capsize=3,label=lab)
ax[0].axhline(0,color='k',lw=1.1,ls='--')
ax[0].set_xlabel(r'$\delta$',fontsize=17); ax[0].set_ylabel(r'$V_\alpha(R=0)$',fontsize=17)
ax[0].set_title('(a)  local contribution',fontsize=17); ax[0].legend(fontsize=12)

# (b) long range, R > 2
for c,lab,col,mk in CH:
    s=h[(h.channel==c)&(h.R>2)].groupby(['delta','seed']).V.mean().reset_index()
    g=s.groupby('delta').V.agg(['mean','sem'])
    ax[1].errorbar(g.index,g['mean'],g['sem'],marker=mk,ms=8,lw=2,color=col,capsize=3,label=lab)
ax[1].axhline(0,color='k',lw=1.1,ls='--')
ax[1].set_xlabel(r'$\delta$',fontsize=17)
ax[1].set_ylabel(r'$\overline{V_\alpha}\ (R>2)$',fontsize=17)
ax[1].set_title('(b)  long-range contribution',fontsize=17); ax[1].legend(fontsize=12)

# (c) V(R) vs R at delta = 0.4 only. R=0 is panel (a) and is 10x larger than the
# rest, so it is excluded here to keep the tail visible on a linear scale.
for c,lab,col,mk in CH:
    s=h[(h.channel==c)&np.isclose(h.delta,0.4)].groupby('R').V.agg(['mean','sem'])
    s=s[(s.index>0)&(s.index<=5)]
    ax[2].errorbar(s.index,s['mean'],s['sem'],marker=mk,ms=7,lw=1.9,color=col,
                   capsize=2.5,label=lab)
ax[2].axhline(0,color='k',lw=1.1,ls='--')
ax[2].axvspan(2.05,5.2,color='0.9',zorder=0)
ax[2].text(3.4,ax[2].get_ylim()[1]*0.82,'long-range window',fontsize=11,
           ha='center',color='0.35')
ax[2].set_xlabel(r'pair separation $R$',fontsize=17)
ax[2].set_ylabel(r'$V_\alpha(R)$',fontsize=17)
ax[2].set_title(r'(c)  $V(R)$ at $\delta=0.4$   ($R=0$ omitted)',fontsize=16)
ax[2].legend(fontsize=12)

for a in ax: a.tick_params(labelsize=13)
fig.suptitle(r'Real-space pairing vertex, $L=12$, $U=4$, half filling',fontsize=19)
plt.tight_layout(rect=[0,0,1,0.93]); plt.savefig('/tmp/prfig.png',dpi=110,facecolor='white')

print('long-range V(R>2), half filling:')
print(f'{"delta":>6}'+''.join(f'{c:>11}' for c,_,_,_ in CH)+'    leader')
for dd in sorted(h.delta.unique()):
    vals={}
    for c,_,_,_ in CH:
        vals[c]=h[(h.channel==c)&(h.R>2)&np.isclose(h.delta,dd)].V.mean()
    lead=max(vals,key=lambda k:vals[k])
    print(f'{dd:6.1f}'+''.join(f'{vals[c]:+11.5f}' for c,_,_,_ in CH)+f'    {lead}')
