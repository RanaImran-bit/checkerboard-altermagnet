#!/usr/bin/env python3
"""Does the d-wave pairing condense at FINITE q under the altermagnet? Track the peak of
P(q) vs tam and t1 (6x6 U=4 doped beta=2, DQMC). tam pushes the FULL-susceptibility peak
off q=0 toward (pi,pi); but the connected VERTEX (true pairing) q=0 piece collapses (and
goes negative), so the finite-q weight is the non-interacting BUBBLE (nesting), not pairing.
t1 keeps the peak at q=0. Honest read: no finite-q *pairing* enhancement here -- the
off-q0 peak is band/nesting. A proper FFLO/PDW test needs the leading pair-vertex
eigenvalue (BSE) per q, not max of connected P(q)."""
import os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
root = os.path.normpath(os.path.join(os.path.dirname(__file__),"..",".."))
lines = [l for l in open(os.path.join(root,"results/dqmc_scan","dqmc_qstar_6x6_b2.txt"))]
blocks = {}
cur=None
for l in lines:
    if l.startswith("#") and "scan" in l: cur = "tam" if "scan tam" in l else "t1"; blocks[cur]=[]
    elif l[0:1].isdigit(): blocks[cur].append([x for x in l.strip().split(",")])
def col(b,i,f=float): return np.array([f(r[i]) for r in blocks[b]])
fig, ax = plt.subplots(1,2,figsize=(11.5,4.6))
for k,(b,ttl) in enumerate((("tam",r"altermagnet $t_{am}$"),("t1",r"anisotropic $t'$ ($t_1$)"))):
    x=col(b,0); fmk=col(b,3); fk0=col(b,4); vmk=col(b,7); vk0=col(b,8)
    a=ax[k]
    a.plot(x,fmk,"C0-o",lw=2,ms=7,label=r"full: peak-$q$ (maxk)")
    a.plot(x,fk0,"C0--s",lw=1.6,ms=6,mfc="none",label=r"full: $q{=}0$")
    a.set_xlabel(x.dtype and f"${'t_{am}' if b=='tam' else 't_1'}$"); a.set_ylabel("full d-wave susceptibility")
    a.set_title(ttl); a.grid(alpha=.3)
    a2=a.twinx()
    a2.plot(x,vmk,"C3-^",lw=2,ms=7,label=r"vertex: peak-$q$")
    a2.plot(x,vk0,"C3:v",lw=1.6,ms=6,label=r"vertex: $q{=}0$")
    a2.axhline(0,color="C3",lw=.6,ls=":"); a2.set_ylabel("vertex (pairing)",color="C3")
    a2.tick_params(axis="y",colors="C3")
    # mark where the peak leaves q=0
    for xi,fm,f0 in zip(x,fmk,fk0):
        if abs(fm-f0)>0.5: a.annotate(r"$q^*\neq 0$",(xi,fm),fontsize=8,color="C0",ha="center",va="bottom")
    h1,l1=a.get_legend_handles_labels(); h2,l2=a2.get_legend_handles_labels()
    a.legend(h1+h2,l1+l2,frameon=False,fontsize=8,loc="lower left")
fig.suptitle(r"Peak momentum of $d$-wave $P(q)$: $t_{am}$ pushes the FULL peak off $q{=}0$ (bubble), "
             r"but the VERTEX collapses (6$\times$6, $\beta{=}2$, $n{\approx}0.85$)",fontsize=10.5)
fig.tight_layout(rect=[0,0,1,0.95])
out=os.path.join(root,"docs","dqmc_qstar_tam_t1.png"); fig.savefig(out,dpi=130); print("wrote",out)
