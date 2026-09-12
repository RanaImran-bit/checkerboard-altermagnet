#!/usr/bin/env python -u
"""Build a real bibliography from the Crossref REST API.

A reference list must not be guessed. Web search summaries give titles and rough
venues but rarely the exact volume and page, and an invented citation is worse
than a missing one. Crossref returns the registered record: authors, container
title, volume, page or article number, year, and DOI, all checkable against the
DOI.

arXiv's API was tried first and rate limits hard (HTTP 429) from this address.

Writes /tmp/refs.json for review, then the curated list is turned into
\\bibitem entries by build_manuscripts.py.
"""
import urllib.request, urllib.parse, json, time, sys

MAIL = "kincmp2022@gmail.com"

QUERIES = [
 # foundations and reviews
 ("Smejkal2022a",  "Beyond conventional ferromagnetism and antiferromagnetism altermagnetism"),
 ("Smejkal2022b",  "Emerging research landscape of altermagnetism"),
 ("SmejkalAHE",    "Anomalous Hall antiferromagnets crystal Hall effect Smejkal"),
 ("Mazin2022",     "Prediction of unconventional magnetism in doped FeSb2 Mazin"),
 ("Bai2024",       "Altermagnetism exploring new frontiers in magnetism and spintronics"),
 ("Jungwirth2024", "Altermagnets and beyond nodal magnetically-ordered phases"),
 # materials
 ("Krempasky2024", "Altermagnetic lifting of Kramers spin degeneracy MnTe"),
 ("Lee2024",       "Broken Kramers degeneracy in altermagnetic MnTe"),
 ("Reimers2024",   "Direct observation of altermagnetic band splitting in CrSb"),
 ("Feng2022",      "An anomalous Hall effect in altermagnetic ruthenium dioxide"),
 ("Fedchenko2024", "Observation of time-reversal symmetry breaking in the band structure of altermagnetic RuO2"),
 ("Jiang2025",     "Metallic altermagnetism KV2Se2O spin resolved ARPES"),
 # theory and interactions
 ("Leeb2024",      "Spontaneous formation of altermagnetism from orbital ordering"),
 ("Duerrnagel2025","Altermagnetic phase transition in a Lieb metal"),
 ("Kaushal2025",   "Altermagnetism in modified Lieb lattice Hubbard model"),
 ("Roig2024",      "Minimal models for altermagnetism"),
 ("Brekke2023",    "Two-dimensional altermagnets superconductivity"),
 ("Ouassou2023",   "dc Josephson effect in altermagnets"),
 ("Zhu2023",       "Topological superconductivity in two-dimensional altermagnetic metals"),
 ("Chakraborty2024","Zero-field finite-momentum Josephson effect altermagnet"),
 # checkerboard and frustrated lattices
 ("Pan2024",       "Superconductivity checkerboard lattice Hubbard model determinant quantum Monte Carlo"),
 ("Runge1992",     "Checkerboard lattice Hubbard model magnetic"),
 ("Fujimoto2002",  "Planar pyrochlore checkerboard lattice antiferromagnet"),
 # methodology
 ("Zhang1997",     "Constrained path Monte Carlo method for fermion ground states"),
 ("ZhangKrakauer2003","Quantum Monte Carlo method using phase-free random walks with Slater determinants"),
 ("BSS1981",       "Monte Carlo calculations of coupled boson-fermion systems"),
 ("Hirsch1985",    "Discrete Hubbard-Stratonovich transformation for fermion lattice models"),
 ("White1989",     "Numerical study of the two-dimensional Hubbard model"),
 ("Loh1990",       "Sign problem in the numerical simulation of many-electron systems"),
 ("Troyer2005",    "Computational complexity and fundamental limitations to fermionic quantum Monte Carlo simulations"),
 ("Qin2016",       "Benchmark study of the two-dimensional Hubbard model with auxiliary-field quantum Monte Carlo"),
 ("LeBlanc2015",   "Solutions of the two-dimensional Hubbard model benchmarks and results from a wide range of numerical algorithms"),
 # Hubbard physics and pairing
 ("Hubbard1963",   "Electron correlations in narrow energy bands"),
 ("Lieb1989",      "Two theorems on the Hubbard model"),
 ("Scalapino2012", "A common thread the pairing interaction for unconventional superconductors"),
 ("Anderson1987",  "The resonating valence bond state in La2CuO4 and superconductivity"),
 ("Arovas2022",    "The Hubbard model annual review"),
 ("Qin2020",       "Absence of superconductivity in the pure two-dimensional Hubbard model"),
 # flat band and Lieb
 ("Lieb1989b",     "Ferrimagnetism flat band Lieb lattice Hubbard"),
 ("Mielke1991",    "Ferromagnetism in the Hubbard model on line graphs"),
 ("Tasaki1992",    "Ferromagnetism in the Hubbard models with degenerate single-electron ground states"),
]

def crossref(q, rows=3):
    url = "https://api.crossref.org/works?" + urllib.parse.urlencode(
        {"query.bibliographic": q, "rows": rows,
         "select": "title,author,container-title,volume,page,issued,DOI,article-number,type"})
    req = urllib.request.Request(url, headers={"User-Agent": f"bibliography/1.0 (mailto:{MAIL})"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)["message"]["items"]

out = {}
for key, q in QUERIES:
    try:
        items = crossref(q)
    except Exception as e:
        print(f"  {key:18s} FAILED {e}", file=sys.stderr); continue
    if not items:
        print(f"  {key:18s} no hit"); continue
    it = items[0]
    au = it.get("author", [])
    names = [f"{(a.get('given','') or '')[:1]}.~{a.get('family','')}" for a in au if a.get("family")]
    rec = dict(key=key,
               authors=names,
               title=(it.get("title") or ["?"])[0],
               journal=(it.get("container-title") or ["?"])[0],
               volume=it.get("volume"),
               page=it.get("page") or it.get("article-number"),
               year=it["issued"]["date-parts"][0][0],
               doi=it["DOI"])
    out[key] = rec
    print(f"  {key:18s} {rec['journal'][:34]:34s} {str(rec['volume']):>5}, {str(rec['page']):>10} ({rec['year']})")
    time.sleep(1.2)

json.dump(out, open("/tmp/refs.json", "w"), indent=1)
print(f"\ncollected {len(out)} of {len(QUERIES)}")
