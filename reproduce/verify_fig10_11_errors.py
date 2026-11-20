#!/usr/bin/env python3
"""Checks the error curves shown in Figs. 10 and 11 OF THIS MANUSCRIPT against the saved responses.
Plotted values: fig10_11_plotted_in_manuscript.csv (line data read from the MATLAB figures that are in images/Results/100x100x60/).
Both methods are compared with the global solver G; phase convention s(phi) = mod(phi + 90, 360) - 180.
Usage: python3 verify_fig10_11_errors.py (needs numpy) -> fig10_11_check.csv"""
import csv, glob, os
import numpy as np
H = os.path.dirname(os.path.abspath(__file__)); RESP = os.path.join(H, "fig10_11_responses"); PD = os.path.join(H, "profile_data")
def rd(p):
    with open(p) as h: skip = int(h.readline().strip().startswith("North"))
    return np.loadtxt(p, skiprows=skip)
def align(d, G, col):
    xs, ys = np.unique(G[:, 0]), np.unique(G[:, 1]); ix = np.abs(d[:, 0, None] - xs).argmin(1); iy = np.abs(d[:, 1, None] - ys).argmin(1)
    look = {(x, y): i for i, (x, y) in enumerate(G[:, :2])}; idx = np.array([look[(xs[a], ys[b])] for a, b in zip(ix, iy)]); return d[np.argsort(idx), col]
S = lambda p: np.mod(p + 90, 360) - 180
e_rho = lambda ref, r: 100 * np.linalg.norm(r - ref) / np.linalg.norm(ref)
e_phi = lambda ref, p: 100 * np.linalg.norm(S(p) - S(ref)) / (np.linalg.norm(np.mod(np.abs(S(ref)), 90)) + 2.2e-16)
cfg = {"001": ("rho_phase_1.00000E-03Hz_100x100x60_40.txt", " 0.0010", "0.10000E-02"), "10": ("rho_phase_1.00000E+01Hz_100x100_60.txt", "10.0000", "0.10000E+02")}
calc = {}
for s, (gn, tok, tk) in cfg.items():
    G = rd(os.path.join(PD, gn))
    for comp in ("xy", "yx"):
        c = 5 if comp == "xy" else 7
        for n in (15, 25, 50, 75, 100, 125, 150, 175, 200):
            g = rd(glob.glob(os.path.join(RESP, f"r{comp}_y_fase_freq={tok}Hz*_np2_NS{n}_tol1.0e-02*"))[0])
            d = rd(glob.glob(os.path.join(RESP, f"GUESSr{comp}_y_fase_freq=*{tk}Hz*_np40_NS{n}_tol1.0e-02*"))[0])
            for meth, resp in (("G-POD", g), ("DD-POD", d)):
                calc[(s, "rho_" + comp, meth, n)] = e_rho(G[:, c], align(resp, G, 2)); calc[(s, "phi_" + comp, meth, n)] = e_phi(G[:, c + 1], align(resp, G, 3))
rows, bad, worst = [], 0, 0.0
for r in csv.DictReader(open(os.path.join(H, "fig10_11_plotted_in_manuscript.csv"))):
    k = (r["suffix"], r["QoI"], r["line"], int(float(r["nPOD"]))); v = calc[k]; pl = float(r["value"]); d = abs(v - pl); worst = max(worst, d); ok = d < 0.0051; bad += not ok
    rows.append([{"001": 0.001, "10": 10.0}[r["suffix"]], r["QoI"], r["line"], k[3], pl, f"{v:.6f}", "OK" if ok else "MISMATCH"])
with open(os.path.join(H, "fig10_11_check.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["frequency_Hz", "QoI", "method(reference G)", "n_POD", "plotted_in_manuscript", "recomputed", "status"]); w.writerows(rows)
print(f"{len(rows) - bad} of {len(rows)} plotted points of Figs. 10/11 reproduce (worst difference {worst:.1e} percentage points).")
