# Cluster-Based Wake Model: Project Overview

This document summarizes the full project for internal orientation before the
repository is uploaded to GitHub. It explains what the OpenFOAM case computes,
what the reduced-order model does, which outputs are important, and which files
are intentionally not part of version control.

## 1. Project Goal

The project studies the wake behind a circular cylinder at

```math
Re_D = \frac{U_\infty D}{\nu} \approx 3900.
```

A three-dimensional LES is generated with OpenFOAM. The simulation data are
then compressed into a low-dimensional, data-driven reduced-order description:

1. run a 3D cylinder wake LES from `t = 0` to `t = 300`,
2. identify the developed unsteady interval from force coefficients,
3. use the window `t = 140 ... 270` for the final ROM analysis,
4. project velocity snapshots onto POD coordinates,
5. cluster the POD trajectory with k-means++,
6. construct CMM/CNM transition matrices,
7. identify characteristic network cycles,
8. compare their frequencies with DMD modes,
9. use deterministic CNM mainly as a limitation diagnostic, and
10. use probabilistic CMM/CNM as a more appropriate branching description.

The repository is meant to contain source code, case setup, and documentation.
Large generated CFD fields, VTK files, NumPy arrays, plots, meshes, and virtual
environments are deliberately ignored by Git.

## 2. OpenFOAM Case

The CFD setup lives in:

```text
openfoam/cases/cylinderLES3D/
```

Important committed files:

```text
0/U
0/p
0/nut
constant/physicalProperties
constant/momentumTransport
system/controlDict
system/fvSchemes
system/fvSolution
system/decomposeParDict
system/createPatchDict
karman_cylinder_3D_structured.geo
Allbuild
Allrun
Allclean
Allkill
Allview
calcRe.py
export_lambda2_pvd.py
lambda2_stats.py
```

The active numerical setup is:

| Quantity | Value |
|---|---:|
| solver | `foamRun` / `incompressibleFluid` |
| turbulence model | LES WALE |
| Reynolds number | `Re_D approx 3900` |
| reference velocity | `U_inf = 1` |
| cylinder diameter | `D = 1` |
| kinematic viscosity | `nu = 1/3900` |
| spanwise length | `L_z = 4D` |
| maximum Courant number | `0.75` |
| maximum time step | `0.01` |
| write interval | `0.25` |
| final time | `300` |

The force coefficient function object writes drag, lift, and moment on the
cylinder patch. The reference area is

```math
A_\mathrm{ref} = D L_z = 4.
```

The drag and lift coefficients are

```math
C_D = \frac{F_x}{\frac{1}{2}\rho U_\infty^2 A_\mathrm{ref}},
\qquad
C_L = \frac{F_y}{\frac{1}{2}\rho U_\infty^2 A_\mathrm{ref}}.
```

## 3. Running the OpenFOAM Case

Load OpenFOAM first:

```bash
source /opt/openfoam13/etc/bashrc
```

Build the mesh and run the case:

```bash
cd openfoam/cases/cylinderLES3D
./Allbuild
./Allrun 4
```

`Allrun 4` updates the decomposition, runs in parallel with four MPI ranks,
and reconstructs the latest time after completion.

To inspect or update Reynolds-number-related parameters:

```bash
./calcRe.py
./calcRe.py --target 3900
```

Use parameter-changing options deliberately, because they modify committed case
configuration files.

## 4. Lambda2 Visualization Data

The Lambda2 criterion identifies vortex cores from the second eigenvalue of

```math
\mathbf S^2 + \boldsymbol{\Omega}^2,
```

where

```math
\mathbf S = \frac{1}{2}(\nabla\boldsymbol u + \nabla\boldsymbol u^T),
\qquad
\boldsymbol\Omega = \frac{1}{2}(\nabla\boldsymbol u - \nabla\boldsymbol u^T).
```

Regions with `lambda2 < 0` are rotation-dominated vortex-core candidates.
For the final visualization window:

```bash
cd openfoam/cases/cylinderLES3D
postProcess -func Lambda2 -time '140:270'
touch cylinderLES3D.foam
paraview cylinderLES3D.foam
```

In ParaView, select `internalMesh`, `U`, and `Lambda2`, move to a developed
time such as `t = 200`, create a `Contour` on `Lambda2`, and start with
`lambda2 = -0.25`. Color the surface by `U` magnitude and adjust camera,
lighting, and background interactively. This is intentionally left to the GUI,
because the view choice is part of making a good explanatory image.

Optional helper scripts:

```bash
pvpython lambda2_stats.py cylinderLES3D.foam
pvpython export_lambda2_pvd.py cylinderLES3D.foam lambda2_series 5 --binary --stats
```

## 5. Final ROM Data Window

The final analysis uses the developed wake interval

```text
t = 140 ... 270
```

This window was chosen from the force coefficients after the full simulation
from `t = 0` to `t = 300`. It captures the main oscillatory wake dynamics while
excluding the early transient and avoiding the entire long-run envelope as a
single modelling target.

Current local input files for the final reference run:

```text
data/processed/forceCoeffs_merged_t0_300.dat
data/processed/pod_U_t140_270.npz
```

Both files are generated artifacts and are ignored by Git.

The final NPZ has the following properties:

| Quantity | Value |
|---|---:|
| time interval | `140 ... 270` |
| sample spacing | `dt = 0.25` |
| samples | `521` |
| stored POD coefficients | `200` |
| schema | v2, with spatial modes and coordinates |

## 6. POD Projection

Let the centered snapshot matrix be

```math
\mathbf X_c =
\begin{bmatrix}
(\boldsymbol x_1 - \bar{\boldsymbol x})^T \\
\vdots \\
(\boldsymbol x_N - \bar{\boldsymbol x})^T
\end{bmatrix}.
```

The VTK-field pipeline uses the method of snapshots: the covariance matrix

```math
\mathbf C = \frac{1}{N-1}\mathbf X_c \mathbf X_c^T
```

is diagonalized first. The temporal POD coefficients are

```math
\boldsymbol a_i(t) = \sigma_i \boldsymbol v_i(t),
```

where `sigma_i` is the singular value and `v_i` is the temporal eigenvector.
The corresponding spatial POD modes are reconstructed blockwise as

```math
\boldsymbol\phi_i = \frac{\mathbf X_c^T \boldsymbol v_i}{\sigma_i}.
```

The final NPZ stores both the coefficient matrix and the spatial mode matrix.
This is why DMD modes can later be projected back into physical space.

## 7. Clustering and Transition Matrices

The POD coefficient trajectory is clustered with k-means++ using

```text
K = 10
```

Cluster labels are reordered for readability, starting with the most populated
cluster and then following frequent outgoing transitions where possible.

For consecutive labels `l_i`, the sampled CMM transition counts are

```math
n_{jk} = \#\{ i : l_i = k,\; l_{i+1} = j \}.
```

The sampled transition matrix is column-stochastic:

```math
P_{jk} = \frac{n_{jk}}{\sum_m n_{mk}}.
```

`P` includes self-transitions and therefore depends on the sampling interval.

For the CNM matrix, repeated labels are compressed into runs. Direct run-to-run
transition counts are written as `n_hat`, and

```math
Q_{jk} = \frac{\widehat n_{jk}}{\sum_{m\neq k}\widehat n_{mk}},
\qquad Q_{kk}=0.
```

The transition-time matrix is based on run residence durations. If run `r`
starts at sample `s_r` and ends at sample `e_r`, then

```math
R_r = t_{e_r} - t_{s_r} + \Delta t_\mathrm{med}.
```

For an observed direct transition `k -> j`,

```math
\tau_r = \frac{R_r + R_{r+1}}{2},
```

and

```math
T_{jk} = \frac{1}{\widehat n_{jk}}\sum_{r:k\to j}\tau_r.
```

Unobserved transitions are stored as `NaN`.

## 8. Deterministic and Probabilistic Dynamics

The deterministic CNM reconstruction follows

```math
j^\star = \arg\max_j Q_{jk}
```

from the current cluster `k` and linearly interpolates between centroids over
`T_{j^\star k}`. This is useful as a diagnostic, but it can collapse onto one
branch in a network with multiple plausible successors.

The probabilistic CMM prediction samples the full column `P[:, k]`, including
self-transitions. It runs many realizations and reports:

```text
predicted label realizations
mean cluster occupation
mean coefficient trajectory
coefficient standard deviation
```

For the final `t = 140 ... 270` run:

| Metric | Value |
|---|---:|
| deterministic occupation L1 error | `0.061420` |
| probabilistic occupation L1 error | `0.028695` |
| number of probabilistic runs | `200` |

The probabilistic occupation is closer to the data, while the deterministic
curve remains useful for showing where a simple single-path CNM is too rigid.

## 9. DMD Comparison

Exact DMD is applied to the same POD coefficient time series. With snapshots as
rows, the implementation forms

```math
\mathbf X_1 = [\boldsymbol x_1,\ldots,\boldsymbol x_{N-1}],
\qquad
\mathbf X_2 = [\boldsymbol x_2,\ldots,\boldsymbol x_N],
```

using transposed arrays internally. A truncated SVD gives

```math
\mathbf X_1 \approx \mathbf U_r \mathbf \Sigma_r \mathbf V_r^*.
```

The reduced linear map is

```math
\widetilde{\mathbf A}
= \mathbf U_r^* \mathbf X_2 \mathbf V_r \mathbf \Sigma_r^{-1}.
```

Exact DMD modes are reconstructed as

```math
\boldsymbol\Phi = \mathbf X_2 \mathbf V_r \mathbf\Sigma_r^{-1}\mathbf W,
```

where `W` contains eigenvectors of `A_tilde`. Discrete eigenvalues `mu` are
converted to continuous rates by

```math
\omega = \frac{\log(\mu)}{\Delta t}.
```

The Strouhal number is

```math
St = \frac{\operatorname{Im}(\omega)}{2\pi}\frac{D}{U_\infty}.
```

For the final run, the CNM cycle frequencies and nearest DMD modes agree well:

| Cycle | CNM reference | DMD mode | Difference |
|---|---:|---:|---:|
| main shedding | `St = 0.203` | `St = 0.203511` | `0.000511` |
| sub-harmonic | `St = 0.136` | `St = 0.136337` | `0.000337` |

This is the strongest result of the project: the sparse transition network
finds physically meaningful cycle frequencies, and DMD independently recovers
matching spectral content.

## 10. Reproducing the Final Result Folder

The current reference folder is:

```text
data/figures/final_t140_270_K10_rerun/
```

Recreate it from the existing local generated files with:

```bash
cd /home/stefan/Uni/cfd_rom_regime_clustering_cylinder
source .venv/bin/activate
export MPLCONFIGDIR=/tmp/matplotlib

OUT="data/figures/final_t140_270_K10_rerun"
NPZ="data/processed/pod_U_t140_270.npz"
FORCE="data/processed/forceCoeffs_merged_t0_300.dat"

mkdir -p "$OUT"

python python/scripts/plot_force_coeffs.py \
    --file "$FORCE" \
    --out "$OUT/00_force_coefficients_t0_300.png" \
    --t-min 0 \
    --t-max 300

python python/scripts/run_cnm_analysis.py \
    --input-type npz \
    --input "$NPZ" \
    --out "$OUT" \
    --K 10 \
    --n-runs 200 \
    --random-state 1

mv "$OUT/summary.json" "$OUT/cnm_summary.json"

python python/scripts/run_dmd_analysis.py \
    --input "$NPZ" \
    --out "$OUT" \
    --dt 0.25

mv "$OUT/summary.json" "$OUT/dmd_summary.json"

python python/scripts/dmd_mode_field.py \
    --npz "$NPZ" \
    --dmd-summary "$OUT/dmd_summary.json" \
    --dmd-dir "$OUT" \
    --mode-st 0.203 \
    --component x \
    --z-slice 2.0 \
    --out "$OUT/14_dmd_mode_field_St0203.png"
```

Important plots:

| File | Use |
|---|---|
| `00_force_coefficients_t0_300.png` | explains transient, overshoot, undershoot, and ROM-window choice |
| `01_pod_energy.png` | documents POD truncation |
| `02_mds_clusters.png` | visualizes clustered state space |
| `04_Q_matrix.png` | main CNM transition probabilities |
| `05_T_matrix.png` | direct transition times |
| `10_dmd_spectrum.png` | DMD/CNM frequency agreement |
| `12_cluster_probabilities_probabilistic.png` | probabilistic improvement over deterministic path |
| `14_dmd_mode_field_St0203.png` | spatial DMD mode for main shedding |

## 11. Interpretation and Limitations

The current model is best described as an interpretable reduced-order analysis,
not as a fully validated predictor.

Key limitations:

- training and reconstruction use the same time interval,
- transition probabilities are empirical and finite-sample estimates,
- deterministic CNM suppresses branching by choosing only `argmax(Q)`,
- probabilistic CMM samples observed transition statistics but does not model
  full physical uncertainty,
- DMD fits one autonomous linear operator over a finite, nonstationary window,
- POD and k-means results depend on the chosen time interval and norm.

Natural next steps:

- second-order CNM, where the next cluster depends on the current and previous
  cluster,
- tCNM or trajectory-based CNM variants with more local dynamical information,
- train/test window validation,
- transition-count diagnostics and uncertainty estimates,
- comparison across `K`, POD truncation, and time-window choices.

## 12. GitHub Upload Checklist

Before committing, check:

```bash
git status --short
git ls-files --others --exclude-standard
```

The following source files should be committed if they are part of the final
project:

```text
python/src/rom/dmd.py
python/scripts/run_dmd_analysis.py
python/scripts/dmd_mode_field.py
openfoam/cases/cylinderLES3D/export_lambda2_pvd.py
openfoam/cases/cylinderLES3D/lambda2_stats.py
```

The following should stay ignored:

```text
.venv/
python/.venv/
data/processed/
data/figures/
openfoam/cases/cylinderLES3D/[0-9]*/
openfoam/cases/cylinderLES3D/constant/polyMesh/
openfoam/cases/cylinderLES3D/postProcessing/
openfoam/cases/cylinderLES3D/VTK/
*.npy
*.npz
*.vtk
*.vtu
*.pvd
*.png
*.msh
```

No generated CFD fields or result arrays are required for GitHub. The project
should be understandable and reproducible from source, with the local generated
artifacts used only for producing figures and checking results.

## 13. References

Arthur, D., and Vassilvitskii, S. (2007). k-means++: The advantages of careful
seeding. *Proceedings of SODA 2007*, 1027-1035.
https://doi.org/10.1145/1283383.1283494

Jeong, J., and Hussain, F. (1995). On the identification of a vortex.
*Journal of Fluid Mechanics*, 285, 69-94.
https://doi.org/10.1017/S0022112095000462

Li, H., Fernex, D., Semaan, R., Tan, J., Morzynski, M., and Noack, B. R.
(2021). Cluster-based network model. *Journal of Fluid Mechanics*, 906, A21.
https://doi.org/10.1017/jfm.2020.785

Nicoud, F., and Ducros, F. (1999). Subgrid-scale stress modelling based on the
square of the velocity gradient tensor. *Flow, Turbulence and Combustion*,
62, 183-200. https://doi.org/10.1023/A:1009995426001

Taira, K., Brunton, S. L., Dawson, S. T. M., Rowley, C. W., Colonius, T.,
McKeon, B. J., Schmidt, O. T., Gordeyev, S., Theofilis, V., and Ukeiley,
L. S. (2017). Modal analysis of fluid flows: An overview. *AIAA Journal*,
55(12), 4013-4041. https://doi.org/10.2514/1.J056060

Tu, J. H., Rowley, C. W., Luchtenburg, D. M., Brunton, S. L., and Kutz,
J. N. (2014). On dynamic mode decomposition: Theory and applications.
*Journal of Computational Dynamics*, 1(2), 391-421.
https://doi.org/10.3934/jcd.2014.1.391
