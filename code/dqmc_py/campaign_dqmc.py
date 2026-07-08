import os, sys
sys.path.insert(0, os.path.expanduser("~/qmc/code/dqmc_py")); os.environ["OMP_NUM_THREADS"]="1"
import numpy as np, multiprocessing as mp
from dqmc import DQMC
# argv: Llist(csv) betalist(csv) OUT
LS=[int(x) for x in sys.argv[1].split(",")]; BS=[float(x) for x in sys.argv[2].split(",")]; OUT=sys.argv[3]
U=4.0; DT=0.0625; NWARM=120; NMEAS=400; NSEED=4; MU=2.0; TAM=0.0  # half-filling, tam=0, sign-free
red=["maxk","k0","r0","rgt"]
def task(a):
    L,beta,seed=a
    q=DQMC(L,L,U,MU,beta,DT,tam=TAM,t1=0.0,seed=seed)
    r=q.run(NWARM,NMEAS,chi=True,kres=True)
    return (L,beta,r["suscV_d"],r["suscV_s"],r["corrV_d"],r["corrV_s"],r["dens"],r["sign"])
if __name__=="__main__":
    mp.set_start_method("fork",force=True)
    jobs=[(L,b,s) for L in LS for b in BS for s in range(1,NSEED+1)]
    with mp.Pool(min(60,len(jobs))) as pool: out=pool.map(task,jobs)
    agg={}
    for L,b,svd,svs,cvd,cvs,dn,sg in out: agg.setdefault((L,b),[]).append((svd,svs,cvd,cvs,dn,sg))
    with open(OUT,"w") as fh:
        fh.write("L,beta,dens,sign,"+",".join(f"suscV_d_{r}" for r in red)+","+",".join(f"suscV_s_{r}" for r in red)
                 +","+",".join(f"corrV_d_{r}" for r in red)+","+",".join(f"corrV_s_{r}" for r in red)+"\n")
        for L in LS:
            for b in BS:
                ds=agg[(L,b)]; dn=np.mean([x[4] for x in ds]); sg=np.mean([x[5] for x in ds])
                vals=[]
                for fi in range(4):
                    for rk in red: vals.append(np.mean([x[fi][rk] for x in ds]))
                fh.write(f"{L},{b:.0f},{dn:.4f},{sg:.4f},"+",".join(f"{v:.4f}" for v in vals)+"\n")
    print("WROTE",OUT)
