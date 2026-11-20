#!/usr/bin/env python3
"""Checks the times plotted in Fig. 8 (total and assembly wall time of G-POD and DD-POD vs n_POD, f = 0.001 and 10 Hz, 100x100x60)
against the simulation logs in timing_logs/.  Plotted values: fig8_plotted_in_manuscript.csv (line data of the MATLAB figures).
Usage: python3 verify_fig8_times.py -> fig8_check.csv.  A value 'matches' when a logged run (TE or TM, any saved series) agrees within 0.011 s;
'typed with fewer digits' means it agrees within 0.025 s (e.g. 7.0 for 7.013)."""
import csv, glob, os, re
H = os.path.dirname(os.path.abspath(__file__)); L = os.path.join(H, "timing_logs")
def parse(path):
    out, cur = [], {}
    for line in open(path, encoding="latin-1"):
        if line.startswith("Mode"): cur = {"mode": line.split()[-1]}
        elif line.startswith("Assembly Time"): cur["asm"] = float(line.split()[-1])
        elif line.startswith("Total Time"): cur["tot"] = float(line.split()[-1]); out.append(dict(cur))
    return out
pat = {"G-POD": "MT_G_ROM/outputs/svdtiming_nprocs{n}_ngx100_ngy100_nsy1_nsz1_T?_freq {t}.txt",
       "DD-POD": "DD_ROM_FREQ/outputs/svdtimings*nprocs{n}_ngx100_ngy100_nsy1_nsz20_T?_freq {t}.txt"}
rows, cnt = [], {"match": 0, "typed with fewer digits": 0, "NO MATCH": 0}
for r in csv.DictReader(open(os.path.join(H, "fig8_plotted_in_manuscript.csv"))):
    f = 0.001 if r["fig"].endswith("001") else 10.0; t = " 0.0010" if f < 1 else "10.0000"; n = int(float(r["x"])); y = float(r["y"])
    meth = "G-POD" if r["line"].startswith("G") else "DD-POD"; key = "asm" if "Assemble" in r["line"] else "tot"
    vals = []
    for p in glob.glob(os.path.join(L, pat[meth].format(n=n, t=t))):
        for b in parse(p): vals.append((b[key], os.path.basename(p)[:36], b["mode"]))
    best = min(vals, key=lambda v: abs(v[0] - y)) if vals else (float("nan"), "", "")
    d = abs(best[0] - y); st = "match" if d <= 0.0111 else ("typed with fewer digits" if d <= 0.0251 else "NO MATCH"); cnt[st] += 1
    rows.append([f, r["line"], n, y, best[0], best[1], best[2], st])
with open(os.path.join(H, "fig8_check.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["frequency_Hz", "line", "n_POD", "plotted_s", "closest_log_s", "log_file", "mode", "status"]); w.writerows(rows)
print(f"{len(rows)} plotted values: {cnt['match']} match a log, {cnt['typed with fewer digits']} match to the digits typed, {cnt['NO MATCH']} do not")
for r in rows:
    if r[-1] == "NO MATCH": print("  NO MATCH:", r[:5])
