#!/usr/bin/env python3
"""Checks the maximum errors quoted in the text of Section 'Real-world model results' (accuracy paragraph) and shows which reference
(G or DD) the plotted DD-POD curves of Figs. 10/11 correspond to.  Usage: python3 verify_text_error_numbers.py (needs numpy)
Text: 'maximum errors across the four QoIs are 3.22% for G-POD and 6.23% for DD-POD at 0.001 Hz (n_POD = 50); 3.06% and 3.69% at 10 Hz;
       0.99% and 4.91% (n_POD = 200, 0.001 Hz); 1.55% and 3.11% (n_POD = 200, 10 Hz)'."""
import glob, os
import numpy as np
RESP, PD = "fig10_11_responses", "profile_data"
def rd(p):
    with open(p) as h: skip = int(h.readline().strip().startswith("North"))
    return np.loadtxt(p, skiprows=skip)
def align(d, G, col):
    xs, ys = np.unique(G[:, 0]), np.unique(G[:, 1]); ix = np.abs(d[:, 0, None] - xs).argmin(1); iy = np.abs(d[:, 1, None] - ys).argmin(1)
    look = {(x, y): i for i, (x, y) in enumerate(G[:, :2])}; idx = np.array([look[(xs[a], ys[b])] for a, b in zip(ix, iy)]); return d[np.argsort(idx), col]
A = lambda p: np.mod(p + 90, 360) - 180
er = lambda ref, r: 100 * np.linalg.norm(r - ref) / np.linalg.norm(ref)
ep = lambda ref, p: 100 * np.linalg.norm(A(p) - A(ref)) / (np.linalg.norm(np.mod(np.abs(A(ref)), 90)) + 2.2e-16)
cfg = {0.001: ("rho_phase_1.00000E-03Hz_100x100x60_40.txt", " 0.0010", "0.10000E-02"), 10.0: ("rho_phase_1.00000E+01Hz_100x100_60.txt", "10.0000", "0.10000E+02")}
text = {(0.001, 50): (3.22, 6.23), (10.0, 50): (3.06, 3.69), (0.001, 200): (0.99, 4.91), (10.0, 200): (1.55, 3.11)}
allok = True
for (f, n), (tg, td) in text.items():
    gname, tok, tk = cfg[f]; G = rd(os.path.join(PD, gname)); out = {"G-POD vs G": [], "DD-POD vs G": [], "DD-POD vs DD": []}
    for comp in ("xy", "yx"):
        c = 5 if comp == "xy" else 7; gr, gp = G[:, c], G[:, c + 1]
        dd = rd(glob.glob(os.path.join(PD, f"r{comp}_*{tk}Hz*_no_svd.txt"))[0]); ddr, ddp = align(dd, G, 2), align(dd, G, 3)
        g = rd(glob.glob(os.path.join(RESP, f"r{comp}_y_fase_freq={tok}Hz*_np2_NS{n}_tol1.0e-02*"))[0])
        d = rd(glob.glob(os.path.join(RESP, f"GUESSr{comp}_y_fase_freq=*{tk}Hz*_np40_NS{n}_tol1.0e-02*"))[0])
        out["G-POD vs G"] += [er(gr, align(g, G, 2)), ep(gp, align(g, G, 3))]
        out["DD-POD vs G"] += [er(gr, align(d, G, 2)), ep(gp, align(d, G, 3))]
        out["DD-POD vs DD"] += [er(ddr, align(d, G, 2)), ep(ddp, align(d, G, 3))]
    m = {k: max(v) for k, v in out.items()}
    ok = abs(m["G-POD vs G"] - tg) < 0.0051 and abs(m["DD-POD vs G"] - td) < 0.0051; allok &= ok
    print(f"f={f:g} Hz, n_POD={n}: text G-POD {tg}, DD-POD {td} | recomputed G-POD vs G {m['G-POD vs G']:.2f}, DD-POD vs G {m['DD-POD vs G']:.2f} (vs DD {m['DD-POD vs DD']:.2f}) -> {'text matches (reference G)' if ok else 'MISMATCH'}")
print("\nAll quoted maxima correspond to errors against G." if allok else "\nSome quoted maxima do not match.")
print("The curves plotted in Figs. 10/11 of this manuscript are the same errors against G (see verify_fig10_11_errors.py: 144 of 144 points reproduce).")
