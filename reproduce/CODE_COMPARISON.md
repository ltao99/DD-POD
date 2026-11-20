# Comparison of the released code with the code that produced the paper's results

**What was done.** The source files of the cluster folders that produced the results of the paper (`MT_DDG_ROM_optimized`: benchmark DD and DD-POD runs and the 100x100x60 DD run at 0.001 Hz; `DD_ROM_FREQ`: 100x100x60 runs at 10 Hz and DD-POD with improved initial guess; `Code_Article`: the article snapshot of December 2025) were compared line by line with `src/` of this release, ignoring comments and white space.
**What was not done.** The released code was **not** re-run and its outputs were **not** compared numerically with the saved responses. The statement below rests on the comparison of the sources only.

| file | MT_DDG_ROM_optimized | DD_ROM_FREQ | Code_Article |
|---|---|---|---|
| `MOD_Constants.f90` | 6 lines differ | 14 lines differ | same |
| `MOD_Quicksort.f90` | same | same | same |
| `calc_error.f90` | 17 lines differ | same | same |
| `commons_3d.f90` | same | 19 lines differ | same |
| `inguess.f90` | same | 4 lines differ | same |
| `l2_utils.f90` | 1 lines differ | same | same |
| `mat_coef_global.f90` | same | same | same |
| `matrix_handling.f90` | same | 36 lines differ | same |
| `max3d-par.f90` | 27 lines differ | 28 lines differ | 10 lines differ |
| `model.f90` | same | 4 lines differ | same |
| `msource.f90` | same | same | same |
| `output.f90` | 3 lines differ | 40 lines differ | same |
| `pri_fields_newz0.f90` | same | same | same |
| `r-b-proc.f90` | 63 lines differ | 22 lines differ | same |
| `read_basis.f90` | 68 lines differ | 116 lines differ | 2 lines differ |
| `readata.f90` | same | 1 lines differ | same |
| `rhs.f90` | 106 lines differ | 11 lines differ | same |
| `set_bet_del_alp.f90` | same | same | same |
| `set_nodes.f90` | 50 lines differ | 69 lines differ | same |
| `set_sub_coordinates.f90` | same | same | same |
| `sizes.f90` | 1 lines differ | 4 lines differ | 3 lines differ |
| `solver_MUMPS.f90` | same | same | same |
| `string_utils.f90` | same | same | same |
| `transm.f90` | 324 lines differ | 12 lines differ | same |

(`sizes.f90` and `MOD_Constants.f90` hold the compile-time settings of each run and differ by design.)

## Nature of the differences

1. **Generalisation to decompositions in y.** The as-run code also allows `nsy > 1`: `set_nodes.f90`, `rhs.f90`, `r-b-proc.f90`, `calc_error.f90` and `transm.f90` contain additional branches with y-direction interfaces (`interf_y`) and block types. For a decomposition without y-division (`nsy = 1`, as in all simulations of the paper) `interf_y` has size 0 (`set_nodes.f90`: `IF (delta_y == 0) ALLOCATE(interf_y(0))`), so these loops are empty. The release contains the `nsy = 1` version. This is an argument from the code, not a numerical test.
2. **Names and instrumentation.** The flag that selects the reduced-order mode is `use_svd` in the as-run code and `use_pod` in the release; output and log file names differ (for example `svdtiming_nprocs<n>_...` with the number of modes in the as-run code, `svdtimings...` in the release); additional timers (for example the mapping time) are written in the as-run code.
3. **`read_basis.f90`.** The as-run versions contain additional options (a module `basis_utils`, handling of empty files, a mean snapshot and lambda files) that were not traced; their effect on the results has not been tested. The file-name fix of release v1.0.1 (`...svd.txt`) is consistent with the names of the basis files of the runs.

## Consequence

The timing logs and responses in `reproduce/` come from the as-run versions. The release is expected to reproduce them for `nsy = 1`, but this has not been verified by a run. A regression test (run the 40x40x36 and 80x80x72 DD and DD-POD cases with this release and compare with `reproduce/benchmark_responses/`) is the way to verify it.
