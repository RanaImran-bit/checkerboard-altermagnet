import os, sys
sys.path.insert(0, os.path.expanduser("~/qmc/code/dqmc_py"))
sys.path.insert(0, os.path.expanduser("~/qmc/code/ftcpqmc_py")); os.environ["OMP_NUM_THREADS"]="1"
import numpy as np, multiprocessing as mp
from ftcpmc import FTCPMC
# argv: L betalist tamlist OUT
L=int(sys.argv[1]); BS=[float(x) for x in sys.argv[2].split(",")]; TS=[float(x) for x in sys.argv[3].split(",")]; OUT=sys.argv[4]
U=4.0; DT=0.0625; NMEAS=180; NW=12; NSEED=3; MU=2.0   # half-filling, constrained (sign=1)
red=["maxk","k0","r0","rgt"]
def task(a):
    beta,tam,seed=a
    q=FTCPMC(L,L,U,MU,beta,DT,tam=tam,t1=0.0,seed=seed,nw=NW,constrained=True,stab=True)
    r=q.run_fb_stab(NMEAS,kres=True,chi=True)
    return (beta,tam,r["suscV_d"],r["suscV_s"],r["corrV_d"],r["corrV_s"],r["dens"],r["sign"])
if __name__=="__main__":
    mp.set_start_method("fork",force=True)
    jobs=[(b,t,s) for b in BS for t in TS for s in range(1,NSEED+1)]
    with mp.Pool(min(30,len(jobs))) as pool: out=pool.map(task,jobs)
    agg={}
    for b,t,svd,svs,cvd,cvs,dn,sg in out: agg.setdefault((b,t),[]).append((svd,svs,cvd,cvs,dn,sg))
    with open(OUT,"w") as fh:
        fh.write("beta,tam,dens,sign,"+",".join(f"suscV_d_{r}" for r in red)+","+",".join(f"suscV_s_{r}" for r in red)
                 +","+",".join(f"corrV_d_{r}" for r in red)+","+",".join(f"corrV_s_{r}" for r in red)+"\n")
        for b in BS:
            for t in TS:
                ds=agg[(b,t)]; dn=np.mean([x[4] for x in ds]); sg=np.mean([x[5] for x in ds])
                vals=[np.mean([x[fi][rk] for x in ds]) for fi in range(4) for rk in red]
                fh.write(f"{b:.0f},{t:.1f},{dn:.4f},{sg:.4f},"+",".join(f"{v:.4f}" for v in vals)+"\n")
    print("WROTE",OUT)
