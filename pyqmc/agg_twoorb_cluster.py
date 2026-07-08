#!/usr/bin/env python3
"""Aggregate the 216-run cluster two-orbital campaign: average over seeds, report
d-wave/ext-s VERTEX correlation (equal-time) and SUSCEPTIBILITY (tau-integrated) vs
anisotropy alpha, for CPQMC (T=0 full model) / DQMC / CP-DQMC (finite-T, U-only), at
L=6x6 and 8x8. Writes tidy per-(method,L) CSVs and a summary figure."""
import os, glob, csv, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SRC = os.path.join(os.path.dirname(__file__), "..", "results", "dqmc_scan", "cluster_twoorb")
OUT = os.path.join(os.path.dirname(__file__), "..", "results", "dqmc_scan")

def load():
    rows = []
    for f in glob.glob(os.path.join(SRC, "*_s*.csv")):
        base = os.path.basename(f)
        meth = base.split("_")[0]
        L = int(base.split("_L")[1].split("_")[0])
        with open(f) as fh:
            for r in csv.DictReader(fh):
                d = {"method": meth, "L": L, "alpha": round(float(r["alpha"]), 3)}
                for k in ("chi_d","chi_d_vertex","chi_s","chi_s_vertex"):
                    d[k] = float(r[k])
                # equal-time vertex correlation: SdV (finite-T) or Cd_tau0/Cs_tau0 (cpqmc full eq-time)
                d["eqd"] = float(r.get("SdV", r.get("Cd_tau0", "nan")))
                d["eqs"] = float(r.get("SsV", r.get("Cs_tau0", "nan")))
                d["sign"] = float(r.get("sign", 1.0))
                rows.append(d)
    return rows

def agg(rows):
    keys = {}
    for r in rows:
        keys.setdefault((r["method"], r["L"], r["alpha"]), []).append(r)
    out = {}
    for (m, L, a), rs in keys.items():
        rec = {"n": len(rs)}
        for col in ("chi_d","chi_d_vertex","chi_s","chi_s_vertex","eqd","eqs","sign"):
            v = np.array([x[col] for x in rs], float); v = v[np.isfinite(v)]
            rec[col] = float(v.mean()) if len(v) else float("nan")
            rec[col+"_se"] = float(v.std(ddof=1)/math.sqrt(len(v))) if len(v) > 1 else 0.0
        out[(m, L, a)] = rec
    return out

def main():
    rows = load(); A = agg(rows)
    methods = ["cpqmc","dqmc","cpdqmc"]; Ls = [6,8]
    titles = {"cpqmc":"CPQMC (T=0, full uxx/uxy/v)","dqmc":"DQMC (finite-T, exact)","cpdqmc":"CP-DQMC (finite-T)"}
    # tidy CSVs + console table
    for m in methods:
        for L in Ls:
            alphas = sorted(a for (mm,LL,a) in A if mm==m and LL==L)
            path = os.path.join(OUT, f"twoorb_{m}_L{L}_agg.csv")
            with open(path,"w",newline="") as f:
                w = csv.writer(f); w.writerow(["alpha","chi_d_vertex","chi_d_vertex_se","chi_s_vertex","chi_s_vertex_se","chi_d","chi_s","eqd","eqs","sign","nseed"])
                for a in alphas:
                    r=A[(m,L,a)]; w.writerow([a,r["chi_d_vertex"],r["chi_d_vertex_se"],r["chi_s_vertex"],r["chi_s_vertex_se"],r["chi_d"],r["chi_s"],r["eqd"],r["eqs"],r["sign"],r["n"]])
    # print headline table: chi_d_vertex (susceptibility) vs alpha
    print("=== d-wave VERTEX SUSCEPTIBILITY chi_dV(alpha), mean +/- se over seeds ===")
    for m in methods:
        for L in Ls:
            alphas = sorted(a for (mm,LL,a) in A if mm==m and LL==L)
            s=" ".join(f"a{a:.1f}:{A[(m,L,a)]['chi_d_vertex']:.2f}±{A[(m,L,a)]['chi_d_vertex_se']:.2f}" for a in alphas)
            print(f"{m:>7} L{L}: {s}")
    # figure: 2 rows (chi vertex susc; equal-time vertex) x 3 cols (methods); 6x6 & 8x8 lines
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    for j,m in enumerate(methods):
        for L,mk,c in [(6,'o','tab:blue'),(8,'s','tab:red')]:
            alphas = sorted(a for (mm,LL,a) in A if mm==m and LL==L)
            if not alphas: continue
            aa=np.array(alphas)
            # row 0: tau-integrated vertex susceptibility (d solid, s dashed)
            yd=[A[(m,L,a)]["chi_d_vertex"] for a in alphas]; yde=[A[(m,L,a)]["chi_d_vertex_se"] for a in alphas]
            ys=[A[(m,L,a)]["chi_s_vertex"] for a in alphas]
            ax[0,j].errorbar(aa,yd,yerr=yde,marker=mk,color=c,label=f"$\\chi_d^V$ {L}x{L}",capsize=2)
            ax[0,j].plot(aa,ys,marker=mk,color=c,ls=':',alpha=0.6,label=f"$\\chi_s^V$ {L}x{L}")
            # row 1: equal-time vertex correlation
            ed=[A[(m,L,a)]["eqd"] for a in alphas]; es=[A[(m,L,a)]["eqs"] for a in alphas]
            ax[1,j].plot(aa,ed,marker=mk,color=c,label=f"eq $d$ {L}x{L}")
            ax[1,j].plot(aa,es,marker=mk,color=c,ls=':',alpha=0.6,label=f"eq $s$ {L}x{L}")
        ax[0,j].set_title(titles[m]); ax[0,j].axhline(0,color='k',lw=0.5)
        ax[0,j].set_xlabel(r"anisotropy $\alpha$"); ax[0,j].set_ylabel(r"vertex susceptibility $\chi^V$")
        ax[1,j].set_xlabel(r"anisotropy $\alpha$"); ax[1,j].set_ylabel("equal-time vertex corr.")
        ax[0,j].legend(fontsize=7); ax[1,j].legend(fontsize=7)
    fig.suptitle("Two-orbital altermagnet: d-wave vs ext-s pairing vs anisotropy (4-seed avg, 6x6 & 8x8)")
    fig.tight_layout()
    p=os.path.join(OUT,"twoorb_cluster_scan.png"); fig.savefig(p,dpi=130)
    print(f"# wrote {p}")

if __name__=="__main__":
    main()
