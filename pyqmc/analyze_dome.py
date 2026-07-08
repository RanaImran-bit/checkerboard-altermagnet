#!/usr/bin/env python3
"""Analyze SC-dome doping results from cluster CSV.

Reads results/cpqmc_nearhalf_L{14,16}.csv (pulled from cluster nodes 250/251)
and concludes whether d-wave leads at the dome doping (n≈0.87).

Usage:
    python pyqmc/analyze_dome.py results/cpqmc_nearhalf_L14.csv results/cpqmc_nearhalf_L16.csv
    python pyqmc/analyze_dome.py results/cpqmc_nearhalf_L14.csv   # single file

Expected CSV columns (from unified_scan --csv):
    tam, suscV_d_{maxk,k0,r0,rgt}, suscV_s_{maxk,k0,r0,rgt},
    corrV_d_{...}, corrV_s_{...}, dens, [pe_lam, pe_dov, pe_sov]

Outputs:
  - Table of d/s ratios per tam
  - Verdict: does d-wave lead at dome doping?
"""
from __future__ import annotations
import sys
import numpy as np

try:
    import csv
except ImportError:
    pass


def load_csv(path):
    rows = []
    with open(path) as fh:
        header = fh.readline().strip().split(",")
        for line in fh:
            vals = line.strip().split(",")
            rows.append({h: v for h, v in zip(header, vals)})
    return header, rows


def to_float(s):
    try:
        return float(s)
    except (ValueError, TypeError):
        return float("nan")


def analyze(path):
    header, rows = load_csv(path)
    print(f"\n=== {path} ===")
    print(f"Rows: {len(rows)}   Columns: {header[:8]}...")

    has_pe = "pe_lam" in header
    has_susc = "suscV_d_k0" in header

    if has_susc:
        print(f"\n{'tam':>6} {'dens':>6} {'suscV_d_k0':>12} {'suscV_s_k0':>12} "
              f"{'d/s ratio':>10} {'verdict':>10}"
              + ("  {'pe_lam':>8} {'pe_dov':>7} {'pe_sov':>7}" if has_pe else ""))
        print("-" * (80 if has_pe else 60))

        tams = sorted(set(to_float(r.get("tam", r.get("TAM", "0"))) for r in rows))
        for tam in tams:
            rr = [r for r in rows if abs(to_float(r.get("tam", "0")) - tam) < 1e-4]
            sd = np.mean([to_float(r["suscV_d_k0"]) for r in rr])
            ss = np.mean([to_float(r["suscV_s_k0"]) for r in rr])
            dn = np.mean([to_float(r.get("dens", "nan")) for r in rr])
            ratio = sd / ss if ss != 0 else float("nan")
            verdict = "d-WAVE" if sd > ss else "ext-s"
            line = (f"{tam:>6.2f} {dn:>6.4f} {sd:>12.5f} {ss:>12.5f} "
                    f"{ratio:>10.3f} {verdict:>10}")
            if has_pe:
                pl = np.mean([to_float(r.get("pe_lam", "nan")) for r in rr])
                pd = np.mean([to_float(r.get("pe_dov", "nan")) for r in rr])
                ps = np.mean([to_float(r.get("pe_sov", "nan")) for r in rr])
                line += f"  {pl:>8.5f} {pd:>7.3f} {ps:>7.3f}"
            print(line)

    print()
    d_wins = sum(1 for r in rows
                 if "suscV_d_k0" in r and
                 to_float(r["suscV_d_k0"]) > to_float(r.get("suscV_s_k0", "0")))
    total = len([r for r in rows if "suscV_d_k0" in r])
    print(f"d-wave leads in {d_wins}/{total} parameter points.")

    if has_pe:
        d_pe = sum(1 for r in rows
                   if to_float(r.get("pe_dov", "0")) > to_float(r.get("pe_sov", "0")))
        print(f"d-wave leads pairing eigenvalue in {d_pe}/{total} points.")


def main():
    if len(sys.argv) < 2:
        print("Usage: python pyqmc/analyze_dome.py results/cpqmc_nearhalf_L14.csv [L16.csv ...]")
        print("  (Pull CSV from cluster first: scp profhokin@thkclusters.duckdns.org:~/qmc/results/cpqmc_nearhalf_L14.csv results/)")
        sys.exit(0)
    for path in sys.argv[1:]:
        try:
            analyze(path)
        except FileNotFoundError:
            print(f"File not found: {path}")
            print("  Pull from cluster: scp -P <port> profhokin@thkclusters.duckdns.org:~/qmc/results/cpqmc_nearhalf_L14.csv results/")


if __name__ == "__main__":
    main()
