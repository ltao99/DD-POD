# Simulation cases of the paper

The cases behind the computing-time and error tables and the figures of Tao and Zyserman (2026) that are verified in `../reproduce/`, with the inputs and the settings they were run with, plus the small 20 x 20 x 18 example shipped with the code. The settings are read from the headers of the run logs and from the names of the output files and agree with the table captions of the paper. The hyperparameter study of the paper (Fig. 5) is not yet documented here.

Inputs: the first three meshes are in `../examples/`; the two large cases are in this folder (the conductivity files are gzip-compressed: `gunzip -k conduc_input.gz`). Checksums are in `SHA256SUMS.txt`.

| Case | Mesh (ngx x ngy x ngz) | Subdomains (nsy x nsz) | MPI processes | Frequency (Hz) | Inputs |
|---|---|---|---|---|---|
| Quick example (not used in the tables of the paper) | 20 x 20 x 18 | see `sizes.f90` of the example | see the example | see `input_mesh_MT` | `../examples/20x20x18/` (+ snapshots in `../data/20x20x18/`) |
| Benchmark DTM1 | 40 x 40 x 36 | 1 x 9 | 18 | 0.001 | `../examples/40x40x36/` |
| Benchmark DTM1 | 80 x 80 x 72 | 1 x 9 | 18 | 0.001 | `../examples/80x80x72/` |
| Benchmark DTM1 | 160 x 160 x 144 | 1 x 36 | 72 | 0.001 | `160x160x144/` |
| Real-world model (Betic Chain, Marti) | 100 x 100 x 60 | 1 x 20 | 40 | 0.001 and 10 | `100x100x60/` |

`nprocs_input = 2 x nsy x nsz` (two processes per subdomain, one per polarization), see the main README.

## Frequency

The frequency is read from `input_mesh_MT` (line 1: number of frequencies, then the frequencies). `100x100x60/input_mesh_MT` is the 0.001 Hz file; `input_mesh_MT_10Hz` is the same file with the frequency line changed to 10 (it is the only difference).

## Parameters of the reduced-order and iterative runs

The names of the output files carry the settings, for example `rxy_y_fase_freq= 0.0010Hz_ngx80_ngy80_ngz72_np18_NS125_tol1.0e-03_relx0.050_penal0.050_maxit999_conv5.0e-04_svd.txt`. The fields correspond to the quantities named in the paper as follows (checked against the captions of Tables 3, 4 and 7: N_s, beta, delta, relaxation, delta_SVD of the benchmark, n_POD and the initial guess of the real-world model):

| Field of the file name | Meaning | Value, benchmark 40 and 80 | Value, 160 x 160 x 144 | Value, real-world model |
|---|---|---|---|---|
| `np` | MPI processes | 18 | 72 | 40 |
| `NS` | number of snapshots (benchmark) / number of modes `n_POD` (real-world model, 15 to 200) | 125 | 125 | 15, 25, 50, 75, 100, 125, 150, 175, 200 |
| `tol` | SVD truncation tolerance `delta_SVD` (value in the real-world runs is read from the file names; the caption of Table 7 does not state it) | 1.0e-3 | 1.0e-3 | 1.0e-2 |
| `relx` | relaxation (0.050 = relaxation at 95 %) | 0.050 | 0.050 | 0.050 |
| `penal` | Robin multiplier, `beta = penal x beta_ref` | 0.050 | 0.050 | 0.050 |
| `maxit` | maximum number of Schwarz iterations | 999 | 999 | 999 |
| `conv` | Schwarz convergence tolerance `delta` | 5.0e-4 | 1.0e-4 | 5.0e-4 |
| prefix `GUESS` | initial guess = converged snapshot nearest in parameter space (paper, real-world model) | - | - | used for DD-POD |

Full-order DD runs (no POD) have the same `np` and `conv`; for the 100 x 100 x 60 model they are the files ending in `_no_svd`, for the 80 x 80 x 72 mesh the DD run of the paper is the file ending in `_beta005` (`beta = 0.05 beta_ref`).

## Which results come from which run

See `../reproduce/INDEX.md`: for every table and figure it names the output file or log that contains the plotted or printed value, and the script that checks it.

## Provenance notes

- The 160 x 160 x 144 conductivity model and mesh are taken from the repository history (commit `498289c`, February 2026); the 100 x 100 x 60 files from the authors' run folder `Results/100x100x60(1airlayer)`. The run folders on the cluster were overwritten after the runs, so byte-identity of these inputs with the files read by the runs of the paper could not be verified by checksum.
- The snapshot matrices and POD bases are not part of the repository (about 30 GB for 80 x 80 x 72 and 200-370 GB per case for the larger meshes). They are produced by full-order DD solves of the snapshot conductivity models; the procedure is described in the main README.
