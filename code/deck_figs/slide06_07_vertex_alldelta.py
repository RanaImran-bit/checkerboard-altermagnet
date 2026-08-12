import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, matplotlib.pyplot as plt

CSV='/home/phd25imran/analysis/data/pairing_master.csv'
OUTDIR='/home/phd25imran/analysis'; L=12; UFIX=4.0; NFIX=1.0
CH=[('chi_son',r'on-site $s$','#7f7f7f','o'),('chi_sext',r'extended $s$','#2ca02c','s'),
    ('chi_d',r'$d_{x^2-y^2}$','#1f77b4','^'),('chi_dxy',r'$d_{xy}$','#d62728','D')]
NC=[c+'_n' for c,_,_,_ in CH]

d=pd.read_csv(CSV); d=d[d.L==L]; d['N']=d.L**2
for c,_,_,_ in CH: d[c+'_n']=d[c]/d.N

def panels(fixcol,fixval,xcol,xlabel,fname,suptitle):
    s=d[np.isclose(d[fixcol],fixval)]
    DS=sorted(s.delta.unique())
    sub={dd: s[np.isclose(s.delta,dd)].groupby(xcol)[NC].agg(['mean','sem']) for dd in DS}
    vals=[g[(c,'mean')] for g in sub.values() for c in NC]
    lo=min(v.min() for v in vals); hi=max(v.max() for v in vals); pad=0.06*(hi-lo)
    ncol=4; nrow=int(np.ceil(len(DS)/ncol))
    fig,ax=plt.subplots(nrow,ncol,figsize=(4.7*ncol,4.2*nrow),facecolor='white',squeeze=False)
    for k,dd in enumerate(DS):
        a=ax[k//ncol][k%ncol]; g=sub[dd]
        for c,lab,col,mk in CH:
            a.errorbar(g.index.values,g[(c+'_n','mean')],g[(c+'_n','sem')],
                       marker=mk,ms=6,lw=1.8,color=col,capsize=2.5,label=lab)
        a.axhline(0,color='k',lw=1.3,ls='--')
        a.set_title(rf'$\delta={dd:g}$',fontsize=16)
        a.set_xlabel(xlabel,fontsize=14); a.tick_params(labelsize=11)
        a.set_ylim(lo-pad,hi+pad)
        if k%ncol==0: a.set_ylabel(r'vertex  $\chi_\alpha/N$',fontsize=14)
    for k in range(len(DS),nrow*ncol): ax[k//ncol][k%ncol].axis('off')
    ax[0][0].legend(fontsize=10,loc='best',framealpha=0.95)
    fig.suptitle(suptitle,fontsize=18)
    plt.tight_layout(rect=[0,0,1,0.94],h_pad=2.5)
    plt.savefig(f'{OUTDIR}/{fname}',dpi=100,facecolor='white'); plt.close()
    return sub

s1=panels('n',NFIX,'U','$U$','G1.png',
   rf'Pairing vertex vs interaction, every anisotropy.  $L={L}$, $n={NFIX:g}$')
s2=panels('U',UFIX,'n','Filling $n$','G2.png',
   rf'Pairing vertex vs filling, every anisotropy.  $L={L}$, $U={UFIX:g}$')

print('FIG 1 (n=1): dxy vertex/N at U=8, per delta:')
print('  '+'  '.join(f'd={k:g}:{v.loc[8.0,("chi_dxy_n","mean")]:+.3f}' for k,v in s1.items()))
print('FIG 2 (U=4): dxy vertex/N at n=1.0 and n=0.778, per delta:')
for tag,nn in [('n=1.000',1.0),('n=0.778',0.7777777777777778)]:
    print(f'  {tag}: '+'  '.join(f'd={k:g}:{v.loc[nn,("chi_dxy_n","mean")]:+.3f}' for k,v in s2.items()))
