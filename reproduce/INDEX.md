# Evidence index: which file proves which displayed result

Status as of 2026-11-20. "Verified" means a number printed in the manuscript was compared, digit by digit, with the raw output of the simulation or with the saved response files, and it agrees. Everything that is **not** verified is listed explicitly in section 6.

Paths starting with `~/old_code/` refer to the cluster (`clonetroop`). The files themselves are copied under `timing_logs/` (same relative paths), so nothing here needs the cluster.

## 1. Timing tables (Tables `perfomance_18procs`, `perfomance_72procs`, `performance_40procs`)

Run `python3 verify_table_values.py`: it reads the logs in `timing_logs/` and writes `table_values_from_logs.csv` (output saved in `verify_table_values_output.txt`).

**17 of 18 printed wall times match a logged run (rounded to 2 decimals).** Both polarizations are stored in the CSV (TE and TM). The one that does not match is in section 6, item 1.

Corrected in this version because they did not match the logs: 40x40x36 G (8.83 -> 8.84), 80x80x72 G (193.78 -> 193.79), 80x80x72 DD (77.01 -> 77.18), 100x100x60 DD-POD at 0.001 Hz (3.16 -> 3.17), and the 40x40x36 G-POD reduction (90 % -> 89 %).

## 2. Timing bar charts (Figs. 6, 7, 9)

Every bar equals the logged time of the step (Assembly, RHS, Analysis/Factorization, Solve, x projection, Mapping, Communication, Total). Source logs:

| Panel | Log (`timing_logs/...`) | Iterations in title |
|---|---|---|
| 6a G 80x80x72 | `MT_G_factorized/outputs/timing_nprocs18_ngx80_ngy80_nsy1_nsz1_TE_freq0.001.txt` (TM block) | - |
| 6b DD | `MT_DDG_ROM_optimized/betas_80/timing_nprocs18_ngx80_ngy80_nsy1_nsz9_TE_freq  0.0010_beta005.txt` | XY 62 (TM file), YX 52 |
| 6c G-POD | `MT_G_ROM/outputs/svdtiming_nprocs2_ngx80_ngy80_nsy1_nsz1_TM_freq  0.0010.txt` | - |
| 6d DD-POD | `MT_DDG_ROM_optimized/betas_80/svdtiming_nprocs18_ngx80_ngy80_nsy1_nsz9_TE_freq  0.0010.txt` | 38 / 38 |
| 7a DD | `MT_DDG_ROM_optimized/betas_160_005/timing_nprocs72_ngx160_ngy160_nsy1_nsz36_TM_freq  0.0010.txt` | XY 232, YX 192 (TE file) |
| 7b DD-POD | `MT_DDG_ROM_optimized/betas_160_001/svdtiming_nprocs72_ngx160_ngy160_nsy1_nsz36_TE_freq  0.0010.txt` | XY 178 (TM file), YX 153 |
| 9a G (10 Hz) | `MT_G_factorized/outputs/timing_nprocs40_ngx100_ngy100_nsy1_nsz1_TE_freq10.txt` (TM block) | - |
| 9b DD (10 Hz) | `DD_ROM_FREQ/outputs/timingsnprocs40_ngx100_ngy100_nsy1_nsz20_TM_freq 10.0000.txt` | XY 147, YX 145 (TE file) |
| 9c G-POD n=200 (10 Hz) | `MT_G_ROM/outputs/svdtiming_nprocs200_ngx100_ngy100_nsy1_nsz1_TM_freq 10.0000.txt` | - |
| 9d DD-POD n=200 (10 Hz) | `DD_ROM_FREQ/outputs/svdtimingsGUESSnprocs200_ngx100_ngy100_nsy1_nsz20_TM_freq 10.0000.txt` | XY 91, YX 92 (TE file) |

Notes
- Fig. 6c: the previous bars (6.18 s, 7.60 s) matched no log; the labels and Table 5 (7.38 s) match the file above, so the bars were rebuilt from it. The log says the G-POD run used 2 MPI processes; the Fig. 6 caption says 18.
- Fig. 9: before this version, panels 9a/9b mixed values of the 0.001 Hz and 10 Hz runs. All four panels now use the 10 Hz runs, as the caption states.
- Fig. 9d previously contained hand-typed numbers (see `figure_scripts/untitled_original_paste_in_script.m`) that matched no log.

## 3. Profiles (Figs. 12 and 13)

Plotted with `figure_scripts/make_profiles.m` from the files in `profile_data/` (profile along x at y = 100 km, 10,000-point response files, all four solvers). For 6 of the 8 panels (Fig. 12 all four, Fig. 13 rho_xy and phi_xy) the curves were compared with the line data stored in the original MATLAB figures: identical (differences below 1e-13). The other two (Fig. 13 rho_yx and phi_yx) have no surviving source figure; they were rebuilt from the same files with the same script and compared with the previous images only visually (not numerically).

| Solver | 0.001 Hz | 10 Hz |
|---|---|---|
| G | `rho_phase_1.00000E-03Hz_100x100x60_40.txt` | `rho_phase_1.00000E+01Hz_100x100_60.txt` |
| G-POD (n=50) | `r?y_y_fase_freq= 0.0010Hz_..._np2_NS50_...svd.txt` | `r?y_y_fase_freq=10.0000Hz_..._np2_NS50_...svd.txt` |
| DD | `r?y_y_fase_freq= 0.10000E-02Hz_..._np40_no_svd.txt` | `r?y_y_fase_freq= 0.10000E+02Hz_..._np40_no_svd.txt` |
| DD-POD (n=50, improved initial guess) | `GUESSr?y_y_fase_freq= 0.10000E-02Hz_..._np40_NS50_...svd.txt` | `GUESSr?y_y_fase_freq= 0.10000E+02Hz_..._np40_NS50_...svd.txt` |

(`?` = x for XY, y for YX.) The y-axis limits are the only thing that changed in this revision.

## 3b. Times versus n_POD (Fig. 8)

Run `python3 verify_fig8_times.py` (uses the 148 G-POD and DD-POD logs in `timing_logs/`; plotted values in `fig8_plotted_in_manuscript.csv`; result in `fig8_check.csv`). **All 72 plotted values trace to a logged run** (69 to within 0.01 s, 3 to the digits typed, e.g. 7.0 for 7.013). Series used: G-POD total and assembly from `MT_G_ROM/outputs/svdtiming_nprocs<n>_..._TM`; DD-POD total from the improved-initial-guess runs (`svdtimingsGUESSnprocs<n>`); **DD-POD assembly from the older runs without the improved guess (`svdtimings_nprocs<n>` / `svdtimingsnprocs<n>`)**, so the DD-POD total and assembly curves come from two different sets of runs (assembly does not depend on the initial guess, so the effect is small).

Corrected in this version (10 Hz panel, `figure_scripts/make_fig8_10Hz.m`): three plotted values had no log and were replaced by the logged ones: G-POD total at n_POD = 175 (14.3 -> 16.11, from `svdtiming_nprocs175_..._TM_freq 10.0000`), G-POD assembly at n_POD = 175 (12.29 -> 14.53, same file) and DD-POD assembly at n_POD = 200 (4.04 -> 3.74, from `svdtimingsnprocs200_..._TE_freq 10.0000`; 4.04 was the 0.001 Hz value). The 'Iterative Scheme' arrow at n_POD = 200 was moved with the marker. Note: the manuscript's 0.001 Hz panel is a later export of the same figure with the 'Iterative Scheme' arrows added; its curves are those of the source figure that was checked.

## 4. Error curves (Figs. 10 and 11, relative L2 errors vs n_POD)

Run `python3 verify_fig10_11_errors.py` (uses `fig10_11_responses/` and `profile_data/`; plotted values read from the figures of this manuscript are in `fig10_11_plotted_in_manuscript.csv`; writes `fig10_11_check.csv`, output in `verify_fig10_11_errors_output.txt`). **All 144 plotted points reproduce** (worst difference 5e-7 percentage points). Both methods are compared with G, as the caption states, and with the phase convention s(phi) = mod(phi + 90, 360) - 180. The figures are drawn by `figure_scripts/make_fig10_11_common_reference_errors.m`. The maxima quoted in the text (3.22/6.23, 3.06/3.69, 0.99/4.91, 1.55/3.11 %) are errors against G and are verified by `verify_text_error_numbers.py`.

History: the Fig. 10/11 images of the previous version (still in the authors' `Code/Results` folder) were different: their DD-POD curves were errors against DD and the G-POD phi_xy curve at 0.001 Hz used another phase wrap. The analysis is kept in `history_previous_fig10_11_images/`; it only matters if a reader compares with the previous version.

## 5. Benchmark error tables (`perfomance_18procs_error` for 40x40x36 and 80x80x72, `perfomance_72procs_error` for 160x160x144)

Run `python3 verify_benchmark_error_table.py` (uses `benchmark_responses/`, writes `benchmark_error_table_check.csv`, output in `verify_benchmark_error_table_output.txt`). **All 28 cells reproduce to two decimals.** The runs behind the rows:

| Mesh | G-POD | DD | DD-POD |
|---|---|---|---|
| 40x40x36 | `np2_NS15_tol1e-03` | `np18` | `np18_NS125_tol1e-03_relx0.050_penal0.050` |
| 80x80x72 | `np2_NS15_tol1e-03` | `np18_beta005` | `np18_NS125_tol1e-03_relx0.050_penal0.050` |
| 160x160x144 (against DD) | - | `np72` | `np72_NS125_tol1e-03_relx0.050_penal0.050_conv1e-04` |

The 80x80x72 DD run (beta = 0.05) is the same run whose timing is in Table 5 and Fig. 6b (`betas_80/..._beta005`); the G-POD NS15 responses are from 2026-01-22, like the G-POD timing logs used in Table 5.

## 6. NOT proven yet (open items)

1. **Table `performance_40procs`, G at 100x100x60, f = 0.001 Hz, 198.13 s: the log is lost, the run is documented; the authors decided on 2026-11-20 not to repeat it.** The timing file of this run was overwritten (the surviving `timing_nprocs40_..._freq0.001.txt`, dated 2026-03-12, reports 299.33 s / 301.33 s and belongs to a later run). Searched without success: all timing logs under `~/old_code`, its moved copy in `/scratch`, and every small text file in `/scratch/<user>`. What exists instead (`slurm_accounting/G_100x100x60_jobs_sacct.txt`): Slurm job 86103 ran from 2025-11-17 00:07:07 to 00:10:27 (200 s elapsed) and wrote the 0.001 Hz G response file at 00:10:26, so 198.13 s of solver time plus about 2 s of start-up is consistent with it, and it is the only G job of the right length. This corroborates but does not reproduce the number. Note that job 86103 was allocated 8 nodes whereas the table caption says two servers. The 10 Hz G time (190.56 s, job 87333) is verified by its log.
2. **Hyperparameter statements** (N_s, beta, delta, relaxation, n_POD values quoted in the text and table captions) were not re-checked against logs. The 'Size LSE' columns of the timing tables were deliberately not checked either (authors' decision, 2026-11-20).
3. **Not checked at all in this pass**: Figs. 1-5 (model, mesh, Pareto of hyperparameters), the 40x40x36 and 80x80x72 profile/error figures other than those listed above, and every number quoted only in the text.
4. The figure scripts (`figure_scripts/*.m`) contain absolute paths of the authors' machine (`/Users/luis/Desktop/...`) and must be given relative paths before they go into the public repository.
