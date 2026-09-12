#!/usr/bin/env python -u
"""Turn the Crossref records into a REVTeX \\bibitem block.

Entries whose Crossref match was wrong are dropped by name rather than silently
kept: a vague bibliographic query sometimes returns a different paper, and a
plausible-looking but incorrect citation is the worst outcome. Two entries are
supplied by hand because they have no usable Crossref record: the group's own
PRB, and a preprint.
"""
import json, sys

R = json.load(open("/tmp/refs.json"))

# Crossref returned the wrong work for these queries. Verified by inspection.
DROP = {"SmejkalAHE", "Jungwirth2024", "Jiang2025", "Pan2024", "Runge1992",
        "Fujimoto2002", "Lieb1989b", "Anderson1987", "Brekke2023x"}
# keys whose record is right but whose name was wrong
RENAME = {"Mazin2022": "Mazin2021", "Hirsch1985": "Hirsch1983",
          "White1989": "Hirsch1985", "Hariki2024": "Hariki2025"}

JOURNAL = {
 "Physical Review X": "Phys. Rev. X", "Physical Review B": "Phys. Rev. B",
 "Physical Review Letters": "Phys. Rev. Lett.", "Physical Review Research": "Phys. Rev. Research",
 "Reviews of Modern Physics": "Rev. Mod. Phys.", "Nature Communications": "Nat. Commun.",
 "Nature Electronics": "Nat. Electron.", "Nature Physics": "Nat. Phys.", "Nature": "Nature",
 "Science": "Science", "Science Advances": "Sci. Adv.",
 "Proceedings of the National Academy of Sciences": "Proc. Natl. Acad. Sci. U.S.A.",
 "Proceedings of the Royal Society of London. Series A. Mathematical and Physical Sciences":
   "Proc. R. Soc. London Ser. A",
 "Annual Review of Condensed Matter Physics": "Annu. Rev. Condens. Matter Phys.",
 "Advanced Functional Materials": "Adv. Funct. Mater.",
 "Europhysics Letters (EPL)": "Europhys. Lett.",
 "Journal of Physics A: Mathematical and General": "J. Phys. A",
}

def fmt(r):
    a = r["authors"]
    if not a: return None
    if len(a) > 4:
        names = a[0] + " \\textit{et al.}"
    elif len(a) == 1:
        names = a[0]
    else:
        names = ", ".join(a[:-1]) + ", and " + a[-1]
    j = JOURNAL.get(r["journal"], r["journal"])
    pg = str(r["page"]).split("-")[0]
    return f"{names},\n{j} \\textbf{{{r['volume']}}}, {pg} ({r['year']})."

lines, kept = [], []
for key, r in sorted(R.items(), key=lambda kv: (kv[1]["year"] or 9999)):
    key = RENAME.get(key, key)
    if key in DROP: continue
    if not (r.get("volume") and r.get("page") and r["journal"] != "?"): continue
    body = fmt(r)
    if not body: continue
    lines.append(f"\\bibitem{{{key}}}\n{body}\n")
    kept.append(key)

# hand entries: no usable Crossref record
lines.append("\\bibitem{OurPRB}\n"
             "Y.~Li \\textit{et al.},\nPhys. Rev. B \\textbf{113}, 134443 (2026).\n")
kept.append("OurPRB")
lines.append("\\bibitem{Pan2024}\n"
             "Y.~Pan, R.~Ma, C.~Chen, Z.~Jia, and T.~Ma,\narXiv:2409.16523.\n")
kept.append("Pan2024")

open("/tmp/bibitems.tex", "w").write("\n".join(lines))
print(f"wrote {len(kept)} bibitems")
print("keys:", " ".join(sorted(kept)))
