import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, glob
import matplotlib.pyplot as plt
TW='/private/tmp/claude-501/-Users-liujiaxin-Desktop-Susceptibility-qmc-platform-master/c2fafe5f-4f29-467d-a412-3fdcefa75ac0/scratchpad/twall'
OUT='/private/tmp/claude-501/-Users-liujiaxin-Desktop-Susceptibility-qmc-platform-master/c2fafe5f-4f29-467d-a412-3fdcefa75ac0/scratchpad/deck2/img'
CH=['chi_son','chi_sext','chi_d','chi_dxy']; N=144
tw=pd.concat([pd.read_csv(f) for f in glob.glob(f'{TW}/*.csv')],ignore_index=True)
tw=tw[tw.nup==72].drop_duplicates(subset=['L','nup','delta','U','seed','apx','apy'])
tw=tw.rename(columns={c+'_vertex':c for c in CH})
per=pd.read_csv('/Users/liujiaxin/Desktop/checkerboard-altermagnet/data/pairing_master.csv')
per=per[(per.L==12)&(per.nup==72)]
def agg(df): return df.groupby(['U','delta'])[CH].agg(['mean','sem'])/N
P,A,B=agg(per),agg(tw[tw.apy==1]),agg(tw[tw.apy==-1])
US=[2.,4.,6.,8.]

fig,ax=plt.subplots(1,4,figsize=(19,4.6),facecolor='white',sharey=True)
for j,u in enumerate(US):
    ds=sorted(set(A.loc[u].index)&set(P.loc[u].index))
    p=[P.loc[(u,x),('chi_dxy','mean')] for x in ds]
    a=[A.loc[(u,x),('chi_dxy','mean')] for x in ds]
    b=[B.loc[(u,x),('chi_dxy','mean')] for x in ds]
    e=[np.sqrt(P.loc[(u,x),('chi_dxy','sem')]**2+(2*A.loc[(u,x),('chi_dxy','sem')])**2
               +B.loc[(u,x),('chi_dxy','sem')]**2)/4 for x in ds]
    avg=[(p[k]+2*a[k]+b[k])/4 for k in range(len(ds))]
    ax[j].plot(ds,p,'--o',ms=7,lw=1.8,color='0.55',label='periodic only')
    ax[j].errorbar(ds,avg,e,fmt='-D',ms=8,lw=2.2,color='#C0392B',capsize=3,
                   label='twist averaged (1:2:1)')
    ax[j].axhline(0,color='k',lw=1.2,ls=':')
    ax[j].set_title(rf'$U={u:g}$',fontsize=18); ax[j].set_xlabel(r'$\delta$',fontsize=17)
    ax[j].tick_params(labelsize=13)
ax[0].set_ylabel(r'$\chi_{d_{xy}}/N$',fontsize=17); ax[0].legend(fontsize=12,loc='lower right')
fig.suptitle(r'$d_{xy}$ pairing vertex: the channel swap survives twist averaging  '
             r'($L=12$, half filling)',fontsize=19)
plt.tight_layout(rect=[0,0,1,0.9]); plt.savefig(f'{OUT}/twist.png',dpi=110,facecolor='white')
print('twist figure written')
for u in US:
    ds=sorted(set(A.loc[u].index)&set(P.loc[u].index))
    r=[(P.loc[(u,x),('chi_dxy','mean')]+2*A.loc[(u,x),('chi_dxy','mean')]+B.loc[(u,x),('chi_dxy','mean')])/4 for x in ds]
    sc=[x for k,x in enumerate(ds) if k and r[k]>0 and r[k-1]<=0]
    print(f'  U={u:g}: crossing at delta = {sc[0] if sc else "n/a"},  peak {max(r):+.4f}')
