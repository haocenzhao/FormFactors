# FormFactors

This repository is used for global fits, uncertainty estimates, result organization, and post-processing plots for nucleon electromagnetic form factors. The current code mainly focuses on the proton and neutron Sachs form factors:

- Proton electric form factor: `G_E^p`
- Proton magnetic form factor: `G_M^p`
- Neutron electric form factor: `G_E^n`
- Neutron magnetic form factor: `G_M^n`

The two main fitting codebases are in `final_fit_proton/` and `final_fit_neutron/`. Final parameterizations, SBS projections, `F_1/F_2`, flavor-separated quark form factors, and related plots are mainly organized in `ff_three_sets/`. The `charge_density/` directory is reserved for later charge-density analysis and is not the main workflow described in this README.

## Directory Structure

```text
FormFactors/
|-- README.md
|-- README_Ch.md
|-- charge_density/
|   `-- charge_distribution.ipynb
|-- final_fit_proton/
|   |-- fitff
|   |-- FF_sumrules_leastsq.py
|   |-- FF_sumrules_coverr.py
|   |-- FF_funcs.py
|   |-- FF_loaddata.py
|   |-- run_coverr.py
|   |-- plot_com.py
|   |-- data/
|   `-- */z*/   # historical outputs for different dates, datasets, and fit options
|-- final_fit_neutron/
|   |-- fitff
|   |-- neut_FF_sumrules_leastsq.py
|   |-- neut_FF_sumrules_coverr.py
|   |-- FF_funcs.py
|   |-- FF_loaddata.py
|   |-- neut_run_coverr.py
|   |-- neut_plot_FF.py
|   |-- data/
|   `-- */z*/   # historical outputs for different dates, datasets, and fit options
`-- ff_three_sets/
    |-- GetFF.py
    |-- GetFF_SBSproj.py
    |-- GetFF_final_SBS.py
    |-- FF_newGMn_SBSproj.py
    |-- FF_newGp_SBSproj.py
    |-- print_gegm_table.py
    |-- plot_new.ipynb
    |-- SBSresult_Apr2026/
    |   |-- proton/
    |   |   `-- SBS_fit_Gp_t0fix7.dat
    |   |-- neutron/
    |   |   `-- out_world_sumrules_leastsq_Q21000_z11_gb5_t0fix7_GEn_SBS.dat
    |   `-- plots/
    |-- plotdata/
    |   |-- GMn_lookup_SBSproj.dat
    |   `-- final_fit_SBS/
    |       |-- SBS_fit_Gp_t0fix7.dat
    |       |-- out_all_sumrules_leastsq_Q21000_z12_gb5_t0fix7_SBS.dat
    |       `-- out_world_sumrules_leastsq_Q21000_z11_gb5_t0fix7_GEn_SBS.dat
    |-- data/
    `-- Quarks_* / Proton_* / Neutron_* plot outputs
```

## Physics and Numerical Method

The core method used in this project is the bounded z expansion. A form factor is written as

```text
G(Q^2) = sum_k a_k z(Q^2)^k
```

with

```text
z = (sqrt(t_cut + Q^2) - sqrt(t_cut - t0))
    / (sqrt(t_cut + Q^2) + sqrt(t_cut - t0))
```

Common settings in the code are:

- `tcut = 4 * mpi**2`
- `t0 = -0.7`, corresponding to `_t0fix7` in file names
- Default proton `kmax = 12`
- Current default neutron entry-point `kmax = 11`
- Coefficient bound usually set to `bound = 5`

The fits include sum-rule constraints to satisfy the normalization at `Q^2 = 0` and the expected large-`Q^2` asymptotic behavior. In the code, `nmax = kmax - 4` is commonly used, meaning that the last several z-expansion coefficients are solved from the sum rules rather than fitted as fully free parameters.

The relation between Sachs form factors and Dirac/Pauli form factors follows the convention in the project notes:

```text
G_E = F_1 - tau * kappa * F_2
G_M = F_1 + kappa * F_2
tau = Q^2 / (4 M_N^2)
```

For later charge-density analysis, note that in a relativistic system one cannot simply interpret `G_E(Q^2)` or `G_M(Q^2)` as the three-dimensional Fourier transform of a static charge or magnetization density. A more robust interpretation usually uses two-dimensional transverse densities in the light-front or infinite-momentum frame.

## Proton Fit: `final_fit_proton/`

`final_fit_proton/` is the main directory for the combined global fit of `G_E^p` and `G_M^p`.

Main files:

- `fitff`: shell wrapper for running proton fits.
- `FF_sumrules_leastsq.py`: main fitting script, using `scipy.optimize.leastsq`.
- `FF_sumrules_coverr.py`: reads the central fit result and generates uncertainty/covariance-propagated output.
- `run_coverr.py`: batch wrapper for `FF_sumrules_coverr.py`.
- `FF_loaddata.py`: loads Mainz, world, polarization, fake high-`Q^2`, SBS, and related data.
- `FF_funcs.py`: basic functions for form factors, cross sections, z expansion, TPE/radiative corrections, and related calculations.
- `plot_com.py` / `plot_com.ipynb`: compares different proton fit versions and makes plots.

### Example

Run from inside `final_fit_proton/`, because many data paths are relative paths:

```bash
cd final_fit_proton
bash fitff 1 1
```

The first `fitff` argument selects the dataset:

```text
1 -> all
2 -> world
3 -> Mainz
```

The second argument selects the fit option:

```text
1 -> default
2 -> no high-Q2 fake constraints, output name contains _noHQ
3 -> no radius constraints, output name contains _noRad
4 -> no extra high-Q2 TPE uncertainty, output name contains _noExtraTPE
5 -> fixed radius option, output name contains _RE8409_RM851
```

Important current defaults in `fitff`:

```text
com    = Apr2026
kmax   = 12
Q2max  = 1000
bound  = 5
mod    = _t0fix7
method = 1
```

The output directory is of the form:

```text
final_fit_proton/Apr2026_all_FF/z12/
```

Fit output files are of the form:

```text
all_sumrules_leastsq_Q21000_z12_gb5_t0fix7.dat
all_sumrules_leastsq_Q21000_z12_gb5_t0fix7.txt
```

Important note: `FF_sumrules_leastsq.py` infers datasets and options from keywords in the output file name. For example, `world`, `all`, and `Mainz` determine `csopt`; `_noHQ`, `_noRad`, and `_noExtraTPE` control high-`Q^2` fake points, radius constraints, and the extra TPE term. Be careful when changing output file names.

### Uncertainty Output

After the central fit is finished, run:

```bash
cd final_fit_proton
python3 run_coverr.py 1 1
```

The first argument selects the dataset, and the second argument selects the output type. `run_coverr.py` generates `out_*` files from the central fit `.dat` files for later post-processing and plotting.

## Neutron Fit: `final_fit_neutron/`

`final_fit_neutron/` is the fitting directory for `G_E^n` and `G_M^n`. Unlike the proton case, the neutron code usually fits `G_E^n` and `G_M^n` separately.

Main files:

- `fitff`: shell wrapper for running neutron fits.
- `neut_FF_sumrules_leastsq.py`: main fitting script.
- `neut_FF_sumrules_coverr.py`: neutron uncertainty/covariance propagation script.
- `neut_run_coverr.py`: batch wrapper for `neut_FF_sumrules_coverr.py`.
- `FF_loaddata.py`: loads neutron world data, fake high-`Q^2` points, and SBS `G_E^n` data.
- `FF_funcs.py`: basic functions used by the neutron fit.
- `neut_plot_FF.py` / `neut_plot_FF.ipynb`: neutron form-factor plotting.

### Example

Run from inside `final_fit_neutron/`:

```bash
cd final_fit_neutron
bash fitff 1   # GEn
bash fitff 2   # GMn
```

The first `fitff` argument is:

```text
1 -> GEn
2 -> GMn
```

Important current defaults in `fitff`:

```text
fitdata = Apr2026_world
kmax    = 11
Q2_Max  = 1000
bound   = 5
fmod    = _t0fix7
method  = 1
```

The output directory is of the form:

```text
final_fit_neutron/Apr2026_world_norm_bound5/z11/
```

Output files are of the form:

```text
Apr2026_world_sumrules_leastsq_Q21000_z11_gb5_t0fix7_GEn.dat
Apr2026_world_sumrules_leastsq_Q21000_z11_gb5_t0fix7_GMn.dat
```

Example neutron uncertainty output:

```bash
cd final_fit_neutron
python3 neut_run_coverr.py 1 1   # GEn default
python3 neut_run_coverr.py 2 1   # GMn default
```

## Summary and Post-Processing: `ff_three_sets/`

`ff_three_sets/` collects the final form-factor parameterizations, SBS projections, `F_1/F_2`, flavor-separated quark form factors, and publication/report-level figures. The current latest Apr2026 SBS fit archive is `SBSresult_Apr2026/`, which should be treated as the primary source for the newest results.

Main files:

- `GetFF.py`: basic form-factor parameterization interface.
- `GetFF_SBSproj.py`: SBS projection parameterization interface, mainly for intermediate results including Tyler's updated `G_M^n`.
- `GetFF_final_SBS.py`: final SBS-updated version, including Haocen's 2026 proton/GEn SBS fit and Tyler's updated `G_M^n`.
- `FF_newGMn_SBSproj.py`: SBS projection interface for Tyler's updated `G_M^n`.
- `FF_newGp_SBSproj.py`: SBS projection interface for updated proton form factors.
- `plot_new.ipynb`: post-processing notebook; the official latest results are the lookup tables and plots under `SBSresult_Apr2026/`.
- `print_gegm_table.py`: calls `GetFF(kID, Q2)` and prints tables for `mu_p G_E^p/G_M^p` and `mu_n G_E^n/G_M^n`.

The `kID` convention for `GetFF(kID, Q2)` is:

```text
1 -> GEp
2 -> GMp
3 -> GEn
4 -> GMn
```

Example:

```bash
cd ff_three_sets
python3 print_gegm_table.py
```

### Apr2026 SBS Results

```text
ff_three_sets/SBSresult_Apr2026/
|-- proton/SBS_fit_Gp_t0fix7.dat
|-- neutron/out_world_sumrules_leastsq_Q21000_z11_gb5_t0fix7_GEn_SBS.dat
`-- plots/
```

The two `.dat` files are lookup tables for the latest Apr2026 fit results:

- `proton/SBS_fit_Gp_t0fix7.dat`: proton FF lookup table, including `Q2`, `z`, `G_E^p`, `G_E^p/G_D`, `G_M^p`, `G_M^p/(mu_p G_D)`, `mu_p G_E^p/G_M^p`, and their uncertainties.
- `neutron/out_world_sumrules_leastsq_Q21000_z11_gb5_t0fix7_GEn_SBS.dat`: neutron `G_E^n` lookup table, including `Q2`, `z`, `G_E^n`, `dG_E^n`, and `G_D`.

The `plots/` directory stores the check plots and summary plots for the Apr2026 fit. Later charge-density / transverse-density analysis should preferentially read the lookup tables in `SBSresult_Apr2026/`, rather than older temporary output directories.

Files named `Quarks_*`, `Proton_*`, and `Neutron_*` in this directory are plot outputs for Sachs form factors, Dirac/Pauli form factors, SBS projections, and up/down quark flavor-separated results.

## Data Flow

A typical workflow is:

```text
data/proton, data/neutron
        |
        v
final_fit_proton/FF_sumrules_leastsq.py
final_fit_neutron/neut_FF_sumrules_leastsq.py
        |
        v
central fit .dat / .txt
        |
        v
FF_sumrules_coverr.py / neut_FF_sumrules_coverr.py
        |
        v
out_* form factor tables and uncertainties
        |
        v
ff_three_sets/GetFF*.py and plot scripts
        |
        v
SBSresult_Apr2026 lookup tables and plots
        |
        v
latest Apr2026 fit products for later density analysis
```

## Dependencies

There is no fixed environment file in this repository. Based on the current scripts, common dependencies include:

- Python 3
- `numpy`
- `scipy`
- `matplotlib`
- Jupyter Notebook, for `.ipynb` files

Some scripts read data files through relative paths, so it is recommended to run scripts from inside the corresponding subdirectory rather than directly from the repository root.

## Important Notes

- Do not casually delete historical output directories. Directories such as `mar30_*`, `may23_*`, `jun14_*`, `oct03_*`, and `Apr2026_*` record results with different dates, datasets, `kmax`, `t0`, constraints, and systematic-uncertainty treatments.
- Do not casually overwrite existing `.dat` files. The main fitting scripts write output in overwrite mode, so rerunning with the same file name will replace old results.
- `ff_three_sets/SBSresult_Apr2026/proton/SBS_fit_Gp_t0fix7.dat` and `ff_three_sets/SBSresult_Apr2026/neutron/out_world_sumrules_leastsq_Q21000_z11_gb5_t0fix7_GEn_SBS.dat` are lookup tables for the latest Apr2026 fit results and will be used in later density analysis.
