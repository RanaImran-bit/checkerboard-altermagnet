#!/usr/bin/env python3
"""Provenance registry for CPQMC runs: which PARAMETERS were computed, with which
CODE VERSION (git commit + driver), at which CLUSTER LOCATION (node + path), and the
resulting observables. SQLite-backed (results/runs.db), stdlib only.

A CAMPAIGN is one scan (a driver + a parameter range, run on a node); a POINT is one
parameter tuple (lx,ly,nup,ndn,U,t0,tam,t1) with its measured observables. One CPMC
run yields many observables, so a point row carries them all (chi_d/chi_d_vtx full &
connected susceptibility, Cd0/Cd0_vtx equal-time, Npp_R0/Npp_Rgt R-resolved, energy).

    python pyqmc/run_db.py init
    python pyqmc/run_db.py add-campaign --name cs_U4 --driver chid_tam_scan.py \
        --node node-258 --path '~/qmc/results/cs' --params 'L8 N31 U4 tam,t1 in 0..0.5' \
        --obs chi_d,chi_d_vtx,Cd0,Cd0_vtx --status running
    python pyqmc/run_db.py ingest --campaign 7 --kind chid file.csv   # or --kind rspace
    python pyqmc/run_db.py list                 # campaigns
    python pyqmc/run_db.py find --U 4 --tam 0.3 --t1 0.0   # is this point done?
    python pyqmc/run_db.py summary -o docs/run_registry.md
"""
from __future__ import annotations
import os, sys, sqlite3, argparse, subprocess, datetime, glob

DB = os.path.join(os.path.dirname(__file__), "..", "results", "runs.db")
OBS_COLS = ["chi_d", "chi_d_vtx", "chi_s", "chi_s_vtx", "Cd0", "Cd0_vtx",
            "Npp_R0", "Npp_Rgt", "energy"]


def _conn():
    os.makedirs(os.path.dirname(DB), exist_ok=True)
    c = sqlite3.connect(DB); c.row_factory = sqlite3.Row; return c


def git_commit():
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"],
                                       cwd=os.path.dirname(__file__)).decode().strip()
    except Exception:
        return "unknown"


def git_commit_of(driver):
    """Commit that last touched the driver script (its code version)."""
    try:
        p = os.path.join(os.path.dirname(__file__), driver)
        return subprocess.check_output(["git", "log", "-1", "--format=%h", "--", p],
                                       cwd=os.path.dirname(__file__)).decode().strip() or git_commit()
    except Exception:
        return git_commit()


def init(_):
    c = _conn()
    c.executescript("""
    CREATE TABLE IF NOT EXISTS campaigns(
      id INTEGER PRIMARY KEY AUTOINCREMENT, name TEXT, driver TEXT, git_commit TEXT,
      node TEXT, result_path TEXT, params TEXT, observables TEXT, status TEXT,
      created_at TEXT, notes TEXT);
    CREATE TABLE IF NOT EXISTS points(
      id INTEGER PRIMARY KEY AUTOINCREMENT, campaign_id INTEGER,
      lx INT, ly INT, nup INT, ndn INT, density REAL,
      U REAL, t0 REAL, tam REAL, t1 REAL, seed TEXT,
      chi_d REAL, chi_d_vtx REAL, chi_s REAL, chi_s_vtx REAL,
      Cd0 REAL, Cd0_vtx REAL, Npp_R0 REAL, Npp_Rgt REAL, energy REAL,
      UNIQUE(campaign_id, lx, ly, nup, ndn, U, t0, tam, t1, seed));
    """)
    c.commit(); print(f"initialized {DB}")


def add_campaign(a):
    c = _conn()
    gc = a.commit or git_commit_of(a.driver)
    cur = c.execute("INSERT INTO campaigns(name,driver,git_commit,node,result_path,params,"
                    "observables,status,created_at,notes) VALUES(?,?,?,?,?,?,?,?,?,?)",
                    (a.name, a.driver, gc, a.node, a.path, a.params, a.obs, a.status,
                     a.date or datetime.date.today().isoformat(), a.notes or ""))
    c.commit(); print(f"campaign {cur.lastrowid}: {a.name}  [{a.driver}@{gc}]  {a.node}:{a.path}")


def _upsert_point(c, cid, d):
    cols = ["campaign_id", "lx", "ly", "nup", "ndn", "density", "U", "t0", "tam", "t1", "seed"] + OBS_COLS
    vals = [cid] + [d.get(k) for k in cols[1:]]
    ph = ",".join("?" * len(cols))
    c.execute(f"INSERT OR REPLACE INTO points({','.join(cols)}) VALUES({ph})", vals)


def ingest(a):
    """Ingest a driver CSV. kind=chid (chid_tam_scan) or rspace (pair_rspace)."""
    c = _conn(); n = 0
    for ln in open(a.file):
        p = ln.split()
        if not (len(p) >= 14 and p[0].isdigit() and p[1].isdigit()):
            continue
        lx, ly, nup, ndn = int(p[0]), int(p[1]), int(p[2]), int(p[3])
        d = dict(lx=lx, ly=ly, nup=nup, ndn=ndn, density=2 * nup / (lx * ly),
                 U=float(p[4]), tam=float(p[5]), t1=float(p[6]), t0=a.t0, seed=p[7])
        if a.kind == "chid":   # lx ly nup ndn U tam t1 seed chi_d e chi_d_vtx e chi_s e chi_s_vtx e Cd0 Cs0 [Cd0_vtx Cs0_vtx]
            d.update(chi_d=float(p[8]), chi_d_vtx=float(p[10]), chi_s=float(p[12]),
                     chi_s_vtx=float(p[14]), Cd0=float(p[16]))
            if len(p) >= 20:
                d["Cd0_vtx"] = float(p[18])
        elif a.kind == "rspace":  # ... seed Npp_R0 e Npp_Rgt e qsum e
            d.update(Npp_R0=float(p[8]), Npp_Rgt=float(p[10]))
        elif a.kind == "unified":  # ...seed (mean err)x: energy Cd0 Cd0_vtx Cs0 Cs0_vtx chi_d chi_d_vtx chi_s chi_s_vtx Npp_R0 Npp_Rgt
            d.update(energy=float(p[8]), Cd0=float(p[10]), Cd0_vtx=float(p[12]),
                     Cs0=float(p[14]), chi_d=float(p[18]), chi_d_vtx=float(p[20]),
                     chi_s=float(p[22]), chi_s_vtx=float(p[24]),
                     Npp_R0=float(p[26]), Npp_Rgt=float(p[28]))
        _upsert_point(c, a.campaign, d); n += 1
    c.commit(); print(f"ingested {n} points into campaign {a.campaign} from {a.file}")


def list_campaigns(_):
    c = _conn()
    rows = c.execute("SELECT c.*, (SELECT COUNT(*) FROM points p WHERE p.campaign_id=c.id) np "
                     "FROM campaigns c ORDER BY id").fetchall()
    print(f"{'id':>3} {'name':<14} {'driver':<18} {'commit':<9} {'node':<10} {'status':<8} {'pts':>5}  params")
    for r in rows:
        print(f"{r['id']:>3} {r['name']:<14} {r['driver']:<18} {r['git_commit']:<9} "
              f"{r['node'] or '-':<10} {r['status']:<8} {r['np']:>5}  {r['params']}")


def find(a):
    c = _conn(); q = "SELECT * FROM points WHERE 1=1"; args = []
    for k in ("lx", "ly", "nup", "U", "tam", "t1"):
        v = getattr(a, k)
        if v is not None:
            q += f" AND abs({k}-?)<1e-6"; args.append(v)
    rows = c.execute(q, args).fetchall()
    if not rows:
        print("no matching points"); return
    print(f"{'L':>3} {'N':>4} {'n':>6} {'U':>4} {'tam':>5} {'t1':>5} | "
          f"{'chi_d_vtx':>10} {'Cd0_vtx':>9} {'Npp_R0':>8}")
    for r in rows:
        print(f"{r['lx']:>3} {r['nup']:>4} {r['density']:>6.3f} {r['U']:>4.1f} {r['tam']:>5.2f} "
              f"{r['t1']:>5.2f} | {_f(r['chi_d_vtx']):>10} {_f(r['Cd0_vtx']):>9} {_f(r['Npp_R0']):>8}")


def _f(x): return f"{x:.4f}" if x is not None else "-"


def summary(a):
    c = _conn(); out = ["# CPQMC run registry\n",
                        f"_auto-generated by pyqmc/run_db.py — {datetime.date.today().isoformat()}_\n"]
    rows = c.execute("SELECT c.*, (SELECT COUNT(*) FROM points p WHERE p.campaign_id=c.id) np "
                     "FROM campaigns c ORDER BY id").fetchall()
    out.append("| id | name | driver | commit | node | path | params | obs | status | points |")
    out.append("|----|------|--------|--------|------|------|--------|-----|--------|--------|")
    for r in rows:
        out.append(f"| {r['id']} | {r['name']} | `{r['driver']}` | `{r['git_commit']}` | "
                   f"{r['node'] or '-'} | `{r['result_path'] or '-'}` | {r['params']} | "
                   f"{r['observables']} | {r['status']} | {r['np']} |")
    txt = "\n".join(out) + "\n"
    if a.out:
        open(a.out, "w").write(txt); print(f"wrote {a.out}")
    else:
        print(txt)


def main():
    ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("init").set_defaults(f=init)
    p = sub.add_parser("add-campaign"); p.set_defaults(f=add_campaign)
    for x in ("name", "driver", "node", "path", "params", "obs", "status", "notes", "commit", "date"):
        p.add_argument(f"--{x}", default=None)
    p = sub.add_parser("ingest"); p.set_defaults(f=ingest)
    p.add_argument("file"); p.add_argument("--campaign", type=int, required=True)
    p.add_argument("--kind", choices=["chid", "rspace", "unified"], required=True)
    p.add_argument("--t0", type=float, default=1.0)
    sub.add_parser("list").set_defaults(f=list_campaigns)
    p = sub.add_parser("find"); p.set_defaults(f=find)
    for x in ("lx", "ly", "nup", "U", "tam", "t1"):
        p.add_argument(f"--{x}", type=float, default=None)
    p = sub.add_parser("summary"); p.set_defaults(f=summary); p.add_argument("-o", "--out", default=None)
    a = ap.parse_args(); a.f(a)


if __name__ == "__main__":
    main()
