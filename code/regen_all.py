"""Re-run every patched figure script and report what regenerated cleanly.

The *_251.py scripts are paste-into-a-notebook copies and end in plt.show(); they
are run here with show() stubbed so the code path is still exercised.
"""
import os, sys, glob, traceback, io, contextlib
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE); sys.path.insert(0, HERE)

SCRIPTS = sorted(f for f in glob.glob("*.py")
                 if ("reggrid" in open(f).read() or "RegularGridInterpolator" in open(f).read())
                 and not f.startswith(("check_", "audit_", "patch_", "regen_", "gridinterp")))

ok, bad = [], []
for f in SCRIPTS:
    plt.close("all")
    src = open(f).read()
    g = {"__name__": "__main__", "__file__": os.path.abspath(f)}
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            plt.show = lambda *a, **k: None
            exec(compile(src, f, "exec"), g)
        ok.append(f)
    except BaseException as e:
        bad.append((f, f"{type(e).__name__}: {e}"))
    finally:
        plt.close("all")

print(f"regenerated cleanly ({len(ok)}/{len(SCRIPTS)}):")
for f in ok: print("   ", f)
if bad:
    print(f"\nfailed ({len(bad)}):")
    for f, e in bad: print(f"    {f}\n        {e}")
