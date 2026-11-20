#!/usr/bin/env python3
"""Checks every wall time printed in the manuscript's timing tables against the raw simulation logs in timing_logs/.
Usage: python3 verify_table_values.py            (writes table_values_from_logs.csv, prints a summary)
A value is 'OK' when it equals the logged Total Time rounded to 2 decimals for the TE or the TM polarization."""
import csv, re, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); L = os.path.join(HERE, "timing_logs")
def totals(path):
    out, mode = {}, None
    for line in open(os.path.join(L, path), encoding="latin-1"):
        if line.startswith("Mode"): mode = line.split()[-1]
        elif line.startswith("Total Time"): out[mode] = float(line.split()[-1])
    return out
def iters(path):
    out, mode = {}, None
    for line in open(os.path.join(L, path), encoding="latin-1"):
        if line.startswith("Mode"): mode = line.split()[-1]
        elif line.startswith("Iterations"): out[mode] = line.split()[-1]
    return out
T = [  # table, mesh, f(Hz), solver, value printed in manuscript, log file
 ("Table perfomance_18procs","40x40x36",0.001,"G",8.84,"MT_G_factorized/outputs/timing_nprocs18_ngx40_ngy40_nsy1_nsz1_TE_freq0.001.txt"),
 ("Table perfomance_18procs","40x40x36",0.001,"G-POD",0.93,"MT_G_ROM/outputs/svdtiming_nprocs2_ngx40_ngy40_nsy1_nsz1_TM_freq  0.0010.txt"),
 ("Table perfomance_18procs","40x40x36",0.001,"DD",2.42,"MT_DDG_ROM_optimized/betas_160_05/timing_nprocs18_ngx40_ngy40_nsy1_nsz9_TM_freq  0.0010.txt"),
 ("Table perfomance_18procs","40x40x36",0.001,"DD-POD",0.30,"MT_DDG_ROM_optimized/betas_160_05/svdtiming_nprocs18_ngx40_ngy40_nsy1_nsz9_TE_freq  0.0010.txt"),
 ("Table perfomance_18procs","80x80x72",0.001,"G",193.79,"MT_G_factorized/outputs/timing_nprocs18_ngx80_ngy80_nsy1_nsz1_TE_freq0.001.txt"),
 ("Table perfomance_18procs","80x80x72",0.001,"G-POD",7.38,"MT_G_ROM/outputs/svdtiming_nprocs2_ngx80_ngy80_nsy1_nsz1_TM_freq  0.0010.txt"),
 ("Table perfomance_18procs","80x80x72",0.001,"DD",77.18,"MT_DDG_ROM_optimized/betas_80/timing_nprocs18_ngx80_ngy80_nsy1_nsz9_TE_freq  0.0010_beta005.txt"),
 ("Table perfomance_18procs","80x80x72",0.001,"DD-POD",1.98,"MT_DDG_ROM_optimized/betas_80/svdtiming_nprocs18_ngx80_ngy80_nsy1_nsz9_TE_freq  0.0010.txt"),
 ("Table perfomance_72procs","160x160x144",0.001,"DD",209,"MT_DDG_ROM_optimized/betas_160_005/timing_nprocs72_ngx160_ngy160_nsy1_nsz36_TM_freq  0.0010.txt"),
 ("Table perfomance_72procs","160x160x144",0.001,"DD-POD",15.74,"MT_DDG_ROM_optimized/betas_160_001/svdtiming_nprocs72_ngx160_ngy160_nsy1_nsz36_TE_freq  0.0010.txt"),
 ("Table performance_40procs","100x100x60",0.001,"G",198.13,"MT_G_factorized/outputs/timing_nprocs40_ngx100_ngy100_nsy1_nsz1_TE_freq0.001.txt"),
 ("Table performance_40procs","100x100x60",0.001,"G-POD n=50",7.93,"MT_G_ROM/outputs/svdtiming_nprocs50_ngx100_ngy100_nsy1_nsz1_TM_freq  0.0010.txt"),
 ("Table performance_40procs","100x100x60",0.001,"DD",23.16,"MT_DDG_ROM_optimized/outputs/timing_nprocs40_ngx100_ngy100_nsy1_nsz20_TM_freq  0.0010.txt"),
 ("Table performance_40procs","100x100x60",0.001,"DD-POD n=50",3.17,"DD_ROM_FREQ/outputs/svdtimingsGUESSnprocs50_ngx100_ngy100_nsy1_nsz20_TM_freq  0.0010.txt"),
 ("Table performance_40procs","100x100x60",10.0,"G",190.56,"MT_G_factorized/outputs/timing_nprocs40_ngx100_ngy100_nsy1_nsz1_TE_freq10.txt"),
 ("Table performance_40procs","100x100x60",10.0,"G-POD n=50",8.06,"MT_G_ROM/outputs/svdtiming_nprocs50_ngx100_ngy100_nsy1_nsz1_TM_freq 10.0000.txt"),
 ("Table performance_40procs","100x100x60",10.0,"DD",24.50,"DD_ROM_FREQ/outputs/timingsnprocs40_ngx100_ngy100_nsy1_nsz20_TM_freq 10.0000.txt"),
 ("Table performance_40procs","100x100x60",10.0,"DD-POD n=50",2.61,"DD_ROM_FREQ/outputs/svdtimingsGUESSnprocs50_ngx100_ngy100_nsy1_nsz20_TM_freq 10.0000.txt"),
]
rows, bad = [], 0
for tab,mesh,f,solver,val,path in T:
    t = totals(path); it = iters(path)
    tol = 0.5 if val == 209 else 0.005 + 1e-9
    hit = [m for m,v in t.items() if abs(v-val) <= tol]
    status = "OK (" + "/".join(hit) + ")" if hit else "MISMATCH"
    if not hit: bad += 1
    rows.append([tab,mesh,f,solver,val,path,t.get("TE",""),t.get("TM",""),it.get("TE",""),it.get("TM",""),status])
with open(os.path.join(HERE,"table_values_from_logs.csv"),"w",newline="") as fh:
    w = csv.writer(fh); w.writerow(["table","mesh","frequency_Hz","solver","manuscript_value_s","log_file(timing_logs/)","log_total_TE_s","log_total_TM_s","iterations_TE","iterations_TM","status"]); w.writerows(rows)
for r in rows: print(f"{r[0][-18:]:18s} {r[1]:12s} f={r[2]:<6} {r[3]:12s} ms={r[4]:<7} TE={r[6]!s:9} TM={r[7]!s:9} {r[10]}")
print(f"\n{len(rows)-bad} of {len(rows)} entries match a logged run; {bad} do not.")
