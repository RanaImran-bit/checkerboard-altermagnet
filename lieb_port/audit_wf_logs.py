# Audit every checkerboard wf.log for the trial-selection bug.
#
# Run() chose the trial by argmin over the SUM OF MEAN-FIELD EIGENVALUES, which
# double counts the interaction and does not rank states by variational energy
# (see lieb_port/patch_hf_selection2.py). This pass does not recompute anything;
# it reports, per run, which seed won and by how much, so we can see
#   (a) runs where a seed other than the physical 'neel' won, and
#   (b) runs where the margin is small enough that the wrong criterion could
#       plausibly have flipped the choice.
#
# Reading only. The checkerboard source is not touched.
import os, re, sys, csv

roots = sys.argv[1:] or [os.path.expanduser("~")]
rows = []
for root in roots:
    for dirpath, dirnames, filenames in os.walk(root):
        if "lieb-altermagnet" in dirpath or "lieb_am" in dirpath:
            dirnames[:] = []
            continue
        if "wf.log" not in filenames:
            continue
        p = os.path.join(dirpath, "wf.log")
        try:
            txt = open(p, errors="replace").read()
        except OSError:
            continue
        att = re.findall(r"attempt\s+(\d+)/(\d+)\s+\(seed=(\w+)\):\s*E\s*=\s*(-?[\d.eE+]+)", txt)
        if not att:
            continue
        pm = re.search(r"Params\s*:\s*(.*)", txt)
        es = [(s, float(e)) for _, _, s, e in att]
        es_sorted = sorted(es, key=lambda x: x[1])
        win, wine = es_sorted[0]
        second, seconde = es_sorted[1] if len(es_sorted) > 1 else ("-", float("nan"))
        neel = [e for s, e in es if s == "neel"]
        rows.append(dict(
            path=p, folder=os.path.basename(dirpath),
            params=(pm.group(1).strip() if pm else ""),
            nseeds=len(es), winner=win, winner_E=wine,
            second=second, second_E=seconde, margin=seconde - wine,
            neel_E=(neel[0] if neel else float("nan")),
            neel_minus_win=((neel[0] - wine) if neel else float("nan")),
        ))

w = csv.DictWriter(sys.stdout, fieldnames=list(rows[0].keys()) if rows else ["path"])
w.writeheader()
for r in rows:
    w.writerow(r)
