# Changelog

## v1.1.0 (2026-11-20)
- Added the inputs of the 160x160x144 benchmark and of the real-world 100x100x60 model (`cases/`, gzip-compressed conductivity files, checksums) and `cases/CASES.md` with the run settings of all cases.
- Added `reproduce/`: output files (surface responses, timing logs, Slurm accounting) from which the tables and figures of the paper are built, scripts that check them, and `INDEX.md` mapping each table and figure to its file.
- README: scope (DD and DD-POD only; G and G-POD results provided as output files), paper reference, reproduction instructions, note on snapshots and bases.
- `Makefile` and the example job script: library and MKL paths are variables (`PATHLIB`, `PATHCOMPILER`, `MKLROOT`, `LIBS`).
- Removed tracked binaries and run logs (`a.out`, `stdout*`, `stderr*`, `time`).
- Added `CITATION.cff` and `reproduce/CODE_COMPARISON.md`.

## v1.0.1 (2026-04-15)
- Fix of the file name of the POD basis files read in `src/read_basis.f90` (`...svd.txt`).

## v1.0.0 (2026-04-09)
- Restructured repository (`src/`, `examples/`, `data/`, `docs/`), documentation and MIT license. Archived on Zenodo on 2026-04-14 (doi:10.5281/zenodo.19567027).
