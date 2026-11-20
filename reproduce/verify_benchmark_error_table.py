#!/usr/bin/env python3
"""Recomputes the relative L2 errors of the benchmark error tables from the saved surface responses in benchmark_responses/
and compares them with the numbers printed in the manuscript (Tables perfomance_18procs_error and perfomance_72procs_error).
Usage: python3 verify_benchmark_error_table.py   (needs numpy)  ->  benchmark_error_table_check.csv

Conventions (identical to the real-world error analysis):
  * every saved surface point is used, matched one-to-one by coordinates (tolerance 1e-3 km), no interpolation;
  * rho error = 100 ||rho - rho_ref|| / ||rho_ref||;
  * phase: s(phi) = mod(phi + 90, 360) - 180; error = 100 ||s(phi) - s(phi_ref)|| / ||mod(|s(phi_ref)|, 90)||.
Columns of the response files: 11 columns for the full-order global solver G (x, y, 1/f, ..., rho_xy, phi_xy, rho_yx, phi_yx),
4 columns (x, y, rho, phi) for every other solver."""
import csv, glob, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); B = os.path.join(HERE, "benchmark_responses")

def read(path):
    with open(path) as h: skip = int(h.readline().strip().startswith("North"))
    return np.loadtxt(path, skiprows=skip)
def vals(d, comp):
    if d.shape[1] >= 9: c = 5 if comp == "xy" else 7; return d[:, [c, c + 1]]
    return d[:, 2:4]
def aligned(d, ref, comp):
    xs, ys = np.unique(ref[:, 0]), np.unique(ref[:, 1]); assert len(xs) * len(ys) == len(ref) == len(d)
    ix = np.abs(d[:, 0, None] - xs).argmin(1); iy = np.abs(d[:, 1, None] - ys).argmin(1)
    assert max(np.abs(d[:, 0] - xs[ix]).max(), np.abs(d[:, 1] - ys[iy]).max()) < 1e-3
    look = {(x, y): i for i, (x, y) in enumerate(ref[:, :2])}
    idx = np.array([look[(xs[a], ys[b])] for a, b in zip(ix, iy)]); assert len(set(idx)) == len(ref)
    return vals(d, comp)[np.argsort(idx)]
def err(resp, ref):
    sh = lambda x: np.mod(x + 90.0, 360.0) - 180.0
    r, p = resp.T; gr, gp = ref.T; dphi = sh(p) - sh(gp)
    return 100 * np.linalg.norm(r - gr) / np.linalg.norm(gr), 100 * np.linalg.norm(dphi) / np.linalg.norm(np.mod(np.abs(sh(gp)), 90.0))
def f1(mesh, pattern):
    hit = glob.glob(os.path.join(B, mesh, pattern)); assert len(hit) == 1, (mesh, pattern, hit); return hit[0]

# mesh, solver, table values (rho_xy, phi_xy, rho_yx, phi_yx), response pattern ({c}=xy|yx), reference
T = [
 ("40x40x36", "G-POD",  [0.47,0.26,0.07,0.03], "r{c}_*np2_NS15_*svd.txt",                    "G"),
 ("40x40x36", "DD",     [0.58,0.35,0.47,0.22], "r{c}_*np18.txt",                             "G"),
 ("40x40x36", "DD-POD", [0.68,0.50,0.54,0.22], "r{c}_*np18_NS125_*penal0.050*svd.txt",       "G"),
 ("80x80x72", "G-POD",  [0.29,0.08,0.06,0.05], "r{c}_*np2_NS15_*svd.txt",                    "G"),
 ("80x80x72", "DD",     [0.53,0.33,0.41,0.20], "r{c}_*np18_beta005.txt",                     "G"),
 ("80x80x72", "DD-POD", [0.65,0.47,0.50,0.23], "r{c}_*np18_NS125_*penal0.050*svd.txt",       "G"),
 ("160x160x144","DD-POD",[0.37,0.17,0.14,0.06],"r{c}_*np72_NS125_*penal0.050*svd.txt",       "DD"),
]
rows, bad = [], 0
for mesh, solver, table, pat, refkind in T:
    out = []; files = []
    for c in ("xy", "yx"):
        if refkind == "G":
            gfile = f1(mesh, "rho_phase_*.txt"); ref_d = read(gfile); ref = vals(ref_d, c); grid = ref_d
        else:
            gfile = f1(mesh, f"r{c}_*np72.txt"); ref_d = read(gfile); ref = vals(ref_d, c); grid = ref_d
        rf = f1(mesh, pat.format(c=c)); files.append(os.path.basename(rf))
        e = err(aligned(read(rf), grid, c), ref); out += list(e)
    ok = [abs(a - b) < 0.0051 for a, b in zip(out, table)]; bad += 4 - sum(ok)
    rows.append([mesh, solver, "G" if refkind == "G" else "DD"] + [f"{x:.4f}" for x in out] + table + ["OK" if all(ok) else "MISMATCH"] + files[:1])
with open(os.path.join(HERE, "benchmark_error_table_check.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["mesh","solver","reference","rho_xy_calc","phi_xy_calc","rho_yx_calc","phi_yx_calc","rho_xy_ms","phi_xy_ms","rho_yx_ms","phi_yx_ms","status","response_file_rxy"]); w.writerows(rows)
for r in rows: print(f"{r[0]:12s} {r[1]:7s} vs {r[2]:2s}  calc {r[3:7]}  ms {r[7:11]}  {r[11]}")
print(f"\n{28-bad} of 28 table cells reproduced to two decimals.")
