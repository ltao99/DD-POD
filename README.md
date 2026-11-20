# Accelerating 3D Magnetotelluric Forward Modelling with Domain Decomposition and Order-Reduction Methods

**Authors/Information:** Luis Tao and Fabio Zyserman  
**Paper reference:** Tao, L. and Zyserman, F.: Accelerating 3D magnetotelluric forward modelling with domain decomposition and order reduction methods, EGUsphere preprint egusphere-2026-700, submitted to Geoscientific Model Development (earlier preprint: https://doi.org/10.31223/X5JV1C).  
**Archive:** every release of this repository is archived on Zenodo; all versions: https://doi.org/10.5281/zenodo.19566290

## Abstract
This repository contains the source code, sample setups, and data for the 3D Magnetotelluric (MT) forward modelling code utilizing Domain Decomposition (DD) and Proper Orthogonal Decomposition (POD) techniques. The purpose of this software is to speed up complex 3D MT simulations by reducing the order of the problem while preserving acceptable accuracy in apparent resistivity and phase. 

## Scope of this release

This repository contains the **domain-decomposition (DD) solver and its reduced-order extension (DD-POD)** that are the subject of the paper. The global full-order solver (G) and the global reduced-order solver (G-POD) that are used in the paper as references are **not** part of this code. Their results (surface responses and timing logs) are provided as output files in `reproduce/`, so that every table and figure can be checked, but they cannot be re-run with this code.

## Technical Context: Domain Decomposition & POD
The computational domain of the 3D MT forward modelling is computationally expensive, forming a large-scale linear system. To efficiently solve it the following methods are employed:

1. **Domain Decomposition (DD):** The code divides the global 3D computational domain into several local subdomains which are solved independently and coupled iteratively at the interfaces. Specifically, this domain decomposition strategy partitions the domain **only in the horizontal directions**, meaning that the subdomains only share horizontal interfaces. This distributes the memory footprint and the computational load across multiple processors (using MPI). One processor is assigned per subdomain and polarization.
2. **Proper Orthogonal Decomposition (POD / Order-Reduction):** In addition to the full-order formulation, the DD solver can utilize Reduced Basis-Proper Orthogonal Decomposition (RB-POD) to accelerate computation. By extracting POD bases from pre-computed high-fidelity snapshots and projecting the solver components into a lower-dimensional space, the scheme drastically decreases the Degrees of Freedom (DoF).

## Inputs

The DD solver requires input files and parameters defining the physical model and numerical setup.

### 1. `sizes.f90` (Model & Solver Configuration)
This file controls the domain partitioning and basic mesh dimensions.
* `nsy, nsz`: Number of subdomains in each direction. (`nsy` is 1).
* `nprocs_input`: Total processors used. Must be `2 × nsy × nsz` (two processors are required per subdomain, one for each polarization).
* `ngx, ngy, ngz`: Number of elements along the X, Y, Z directions. `ngz` must be a multiple of `nsz`.

### 2. `input_mesh_MT` (Discretization & Frequencies)
Defines the spatial discretization and frequency sampling. The file structure starts with the number of frequencies and their values. Afterward, for each dimension (X, Y, Z), it defines coordinates of the macro-nodes in meters and the number of partitions between them.

### 3. `conduc_input` (Conductivity Model)
Defines the conductivity distribution (in S/m) of each element in the 3D mesh (X × Y × Z dimensions). The file appends a 1D vertically-layered background model, consisting of layer thicknesses (where the last layer has no thickness assigned) and layer conductivities. Air layers must be explicitly included.

### 4. `MOD_Constants.f90` (Physical Constants & Hyperparameters)
Controls full-order or reduced-order flags and DD hyper-parameters:
* `use_pod`: Set to `.TRUE.` for Reduced-Order Modelling, `.FALSE.` for Full-Order.
* **POD settings**: `Nsnapshots` (number of snapshots), `tol_svd` (SVD truncation tolerance), `pod_path` (directory for snapshot matrices folder.
* **DD settings**: `relx` (relaxation parameter), `penal` (penalization factor), `maxiter`, and `convergence_tol`.

## Outputs

The solver evaluates the physical quantities and dumps them to the outputs directory:
1. **Electromagnetic Fields:** Electric and magnetic field components evaluated at the midpoints of the surface elements stored in the /outputs.
2. **Apparent Resistivity & Phase:** Computed for each polarization evaluated at the midpoints of the surface mesh stored in /outputs.
3. **Performance Logs:** Runtime metrics, timings for each computational phase, and iterative convergence progress stored in /outputs.

## Repository Structure & Data
* **`src/`**: Fortran source code `*.f90`.
* **`examples/`**: Contains configurations (`conduc_input`, `input_mesh_MT`), scripts (`job_mt_monolithic.sh`), and standard parameter setups for 20x20x18, 40x40x36, and 80x80x72 mesh executions. Those example folders define the inputs to generate full-order models.
* **`data/`**: Stores the snapshots databases. `data/20x20x18/Snapshots` contains the high-fidelity state snapshots required to deploy the POD-based reduction.
* **`cases/`**: inputs (conductivity model, mesh, frequency) and the run settings of the two largest cases of the paper (160x160x144 and the real-world 100x100x60 model), and `CASES.md`, which lists the settings of all cases. The conductivity files are gzip-compressed.
* **`reproduce/`**: the small output files (surface responses, timing logs, Slurm accounting) from which the tables and figures of the paper are built, the scripts that check them, and `INDEX.md`, which names for every table and figure the file that contains its values.
* **`docs/`**: Includes the full PDF `Code_Guide.pdf` for extended context.

**Snapshots and POD bases.** Only the snapshots of the 20x20x18 example are included. The snapshot matrices and bases of the larger cases are far too large to archive (about 30 GB for 80x80x72 and 200-370 GB per case for 100x100x60 and 160x160x144). They are produced by full-order DD solves (`usepod = .FALSE.`) of the snapshot conductivity models; the files that are read by the reduced-order runs are named `sol_<TE|TM>_<subdomain>svd.txt` (see `src/read_basis.f90`).

## Note on the code version

The results of the paper were produced with development versions of this code that also allowed a decomposition along y (`nsy > 1`) and used other output file names. For the decomposition used in the paper (`nsy = 1`) the additional y-interface loops are empty. The comparison of the sources is documented in `reproduce/CODE_COMPARISON.md`; the release has not been re-run against the saved responses.

## Reproducing the tables and figures of the paper

```bash
cd reproduce
python3 verify_table_values.py            # computing-time tables against the run logs
python3 verify_benchmark_error_table.py    # benchmark error tables against the saved responses
python3 verify_fig8_times.py               # Fig. 8
python3 verify_fig10_11_errors.py          # Figs. 10 and 11
python3 verify_text_error_numbers.py       # error maxima quoted in the text
```

These scripts need Python 3 with `numpy` and read only the files in `reproduce/`. `reproduce/INDEX.md` lists, for each table and figure, the file that contains the plotted or printed value, and states what is **not** covered: one computing time (G at 100x100x60 and 0.001 Hz) has no log file, and Figs. 1-5 and the hyperparameter study are not covered. The MATLAB scripts in `reproduce/figure_scripts/` draw the figures; they still contain absolute paths of the authors' machine.

## Prerequisites
- **Fortran Compiler:** Intel `mpif90` / Intel oneAPI
- **MPI:** For distributed Domain Decomposition
- **MKL (Math Kernel Library):** For BLAS and LAPACK routines
- **MUMPS & PARMETIS/METIS:** For the distributed sparse direct solver

## Installation & Usage
1. Configure build environment: the `Makefile` paths are variables with defaults; set them on the command line or in the environment, e.g. `make PATHLIB=/path/to/libs MKLROOT=/path/to/mkl` (`PATHLIB` must contain the MUMPS, METIS/ParMETIS and ScaLAPACK libraries; MUMPS 5.1.2 was used for the paper).
2. Build the program: `make`
3. Execute via MPI. For instance, `nsy=1` and `nsz=3` uses 6 processes:
```bash
mpirun -np 6 ./max3d-par
```
*(Ensure that all required input files the code targets, such as `conduc_input`, `input_mesh_MT`, and others are copied to the working directory alongside the executable).*

## Acknowledgments
Development was carried out as part of the Marie Skłodowska-Curie EarthSafe project (EU Horizon Europe programme grant No. 101120556). Special thanks to Prof. Fabio Zyserman for the foundation code and Dr. Constanza Manassero for stiffness matrix assemblage routines.
