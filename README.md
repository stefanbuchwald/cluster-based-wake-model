# CFD wake reduced-order model

A small research and teaching project for building a cluster-based
reduced-order model (ROM) of a cylinder wake.

The repository contains source code and configuration only. Simulation fields,
VTK files, POD matrices, plots, and other generated results are deliberately
excluded from Git.

The workflow is:

```text
OpenFOAM snapshots or force coefficients
                  |
                  v
        centering and POD/SVD
                  |
                  v
           k-means++ clusters
                  |
                  v
     transition matrices P, Q, and T
                  |
                  v
       deterministic CNM trajectory
```

The OpenFOAM case is documented separately in
[`openfoam/README.md`](openfoam/README.md). This document concentrates on the
data analysis and its mathematics.

## Repository layout

```text
.
|-- openfoam/
|   |-- README.md
|   `-- cases/cylinderLES3D/
|-- python/
|   |-- scripts/
|   |   |-- field_pod_from_vtk.py
|   |   |-- plot_force_coeffs.py
|   |   `-- run_cnm_analysis.py
|   `-- requirements.txt
|-- scripts/sweep_cnm_windows.sh
|-- postprocessing.sh
`-- package_project.sh
```

The files below are generated locally and are not part of the repository:

- OpenFOAM time directories, `processor*`, `postProcessing/`, and `VTK/`
- Gmsh meshes and OpenFOAM `polyMesh`
- `data/processed/`, `data/figures/`, and `data/results/`
- NumPy arrays, plots, videos, logs, virtual environments, and ZIP archives

## Installation

Create the Python environment in the project root:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r python/requirements.txt
```

The force-coefficient analysis needs NumPy, pandas, Matplotlib, SciPy, and
scikit-learn. Reading OpenFOAM VTK fields additionally needs PyVista and VTK.

## Inputs

### Force coefficients

`run_cnm_analysis.py --input-type force` reads an OpenFOAM
`forceCoeffs.dat` or `coefficient.dat` file. By default, the state vector is

$$
\boldsymbol{x}(t_i)
=
\begin{bmatrix}
C_D(t_i) & C_L(t_i) & C_M(t_i)
\end{bmatrix}.
$$

Because these observables have different scales, every column is standardized:

$$
\widetilde{x}_{ij}
=
\frac{x_{ij}-\mu_j}{s_j},
$$

where $\mu_j$ and $s_j$ are the sample mean and standard deviation of feature
$j$.

### Velocity or pressure fields

`field_pod_from_vtk.py` reads cell data from a sequence of VTK snapshots.
Velocity components and optional pressure values are flattened into one
snapshot vector:

$$
\boldsymbol{x}_i
=
\begin{bmatrix}
U_x & U_y & U_z & p
\end{bmatrix}_{t_i}^{\mathsf T}.
$$

The pressure block can be multiplied by `--p-scale`. This matters because a
combined velocity-pressure norm otherwise mixes quantities with different
units and physical meanings.

The complete snapshot matrix is stored locally as a memory-mapped `.npy` file.
It is intentionally ignored by Git because it can be several gigabytes large.

## Mathematical method

Let $N$ be the number of snapshots and $M$ the number of input features. The
data matrix uses one snapshot per row:

$$
\mathbf{X}
=
\begin{bmatrix}
\boldsymbol{x}_1^{\mathsf T}\\
\vdots\\
\boldsymbol{x}_N^{\mathsf T}
\end{bmatrix}
\in \mathbb{R}^{N\times M}.
$$

### 1. Mean subtraction and POD

First, the temporal mean is removed:

$$
\overline{\boldsymbol{x}}
=
\frac{1}{N}\sum_{i=1}^{N}\boldsymbol{x}_i,
\qquad
\mathbf{X}'_{i,:}
=
\boldsymbol{x}_i^{\mathsf T}
-\overline{\boldsymbol{x}}^{\mathsf T}.
$$

For a small feature matrix, `run_cnm_analysis.py` computes the singular value
decomposition directly:

$$
\mathbf{X}'
=
\mathbf{U}\boldsymbol{\Sigma}\mathbf{V}^{\mathsf T}.
$$

The columns of $\mathbf{V}$ are spatial or feature-space POD modes. The
temporal POD coefficients used for clustering are

$$
\mathbf{A}
=
\mathbf{U}_r\boldsymbol{\Sigma}_r
\in \mathbb{R}^{N\times r}.
$$

For large VTK fields, `field_pod_from_vtk.py` uses the method of snapshots. It
forms the much smaller temporal covariance matrix

$$
\mathbf{C}
=
\frac{1}{N-1}\mathbf{X}'\mathbf{X}'^{\mathsf T}
\in \mathbb{R}^{N\times N}
$$

and solves

$$
\mathbf{C}\boldsymbol{\psi}_k
=
\lambda_k\boldsymbol{\psi}_k.
$$

The corresponding singular values and temporal coefficients are

$$
\sigma_k=\sqrt{(N-1)\lambda_k},
\qquad
a_k(t_i)=\sigma_k\psi_{ik}.
$$

The modal energy fractions are

$$
e_k=\frac{\sigma_k^2}{\sum_\ell \sigma_\ell^2},
\qquad
E_r=\sum_{k=1}^{r}e_k.
$$

The smallest $r$ satisfying the requested energy threshold is selected, unless
`--n-modes-max` imposes a lower hard limit. Therefore, requesting
`--energy-threshold 0.9999 --n-modes-max 30` does **not** guarantee that 30
modes contain 99.99% of the energy. In the existing 161-snapshot field run,
the first 30 velocity modes contained about 89.27%.

When a POD `.npz` file is passed to `run_cnm_analysis.py`, the current
implementation centers its coefficient matrix and performs another SVD before
clustering. If all supplied dimensions are retained, this is an orthogonal
change of coordinates and preserves Euclidean distances. An additional mode
limit would perform a second truncation.

### 2. k-means++ clustering

Each row $\boldsymbol{a}_i$ of $\mathbf{A}$ is assigned to one of $K$ clusters.
The centroids $\boldsymbol{c}_k$ minimize the within-cluster variance

$$
J
=
\sum_{k=0}^{K-1}
\sum_{i:\,l_i=k}
\left\|
\boldsymbol{a}_i-\boldsymbol{c}_k
\right\|_2^2.
$$

The implementation uses k-means++ initialization, 50 restarts, and a fixed
random seed. After fitting, cluster numbers are reordered for readability:
the most populated cluster becomes $C_0$, followed where possible by its most
frequent outgoing transition. This changes labels, not the clustering.

The empirical occupation probability is

$$
p_k
=
\frac{1}{N}\sum_{i=1}^{N}\mathbb{1}(l_i=k).
$$

### 3. CMM transition matrix

For every pair of consecutive samples, define

$$
n_{jk}
=
\#\left\{i:l_i=k,\ l_{i+1}=j\right\}.
$$

The cluster-based Markov model (CMM) matrix is column-stochastic:

$$
P_{jk}
=
\frac{n_{jk}}{\sum_m n_{mk}}.
$$

Thus, column $k$ contains the probability of the next sampled state, given
that the current sampled state is $C_k$. Its diagonal contains persistence
within a cluster.

### 4. Direct-transition matrix

For the cluster-based network model, repeated labels are compressed into runs:

```text
0, 0, 0, 2, 2, 1, 1  ->  0, 2, 1
```

Let $\widehat n_{jk}$ count transitions between consecutive runs. Self
transitions no longer occur. The direct-transition matrix is

$$
Q_{jk}
=
\frac{\widehat n_{jk}}
{\sum_{m\ne k}\widehat n_{mk}},
\qquad Q_{kk}=0.
$$

Columns with no observed departure remain zero.

### 5. Residence and transition times

Let run $r$ start at sample $s_r$ and end at sample $e_r$. With

$$
\Delta t_{\mathrm{med}}
=
\operatorname{median}(t_{i+1}-t_i),
$$

the implemented residence duration is

$$
R_r
=
t_{e_r}-t_{s_r}+\Delta t_{\mathrm{med}}.
$$

For an observed run transition $k\rightarrow j$, the individual transition
time is

$$
\tau_r
=
\frac{R_r+R_{r+1}}{2}.
$$

The matrix entry is the conditional average over all observed
$k\rightarrow j$ transitions:

$$
T_{jk}
=
\frac{1}{\widehat n_{jk}}
\sum_{r:\,k\rightarrow j}\tau_r.
$$

Unobserved transitions are stored as `NaN`.

### 6. Deterministic CNM reconstruction

Starting from cluster $k$, the current implementation chooses the most likely
direct successor:

$$
j^\star
=
\operatorname*{arg\,max}_{j}Q_{jk}.
$$

During the interval

$$
t_{\mathrm{arrive}}-t_{\mathrm{depart}}
=
T_{j^\star k},
$$

the state is linearly interpolated between centroids:

$$
\widehat{\boldsymbol{a}}(t)
=
(1-\alpha)\boldsymbol{c}_k
+\alpha\boldsymbol{c}_{j^\star},
\qquad
\alpha
=
\frac{t-t_{\mathrm{depart}}}{T_{j^\star k}}.
$$

The displayed cluster label switches from $k$ to $j^\star$ at
$\alpha=0.5$. If a cluster has no observed successor, the predicted state
remains at its centroid.

This reconstruction is deterministic. It does not sample from $Q$ and
therefore cannot reproduce stochastic branching or provide uncertainty
estimates.

### 7. MDS visualization

Multidimensional scaling embeds snapshots and centroids into a three-dimensional
plot by approximately preserving pairwise Euclidean distances:

$$
\min_{\boldsymbol{y}_i}
\sum_{i<j}
\left(
\|\boldsymbol{y}_i-\boldsymbol{y}_j\|_2
-d_{ij}
\right)^2.
$$

MDS is used only for visualization. It does not affect clustering or the
transition matrices.

## Running the analysis

From force coefficients:

```bash
python python/scripts/run_cnm_analysis.py \
    --input-type force \
    --input openfoam/cases/cylinderLES3D/postProcessing/forceCoeffs/0/forceCoeffs.dat \
    --out data/figures/cnm_force_K10 \
    --K 10 \
    --t-min 60
```

From field POD coefficients:

```bash
python python/scripts/field_pod_from_vtk.py \
    --vtk-glob "openfoam/cases/cylinderLES3D/VTK/cylinderLES3D_*.vtk" \
    --fields U \
    --out data/processed/pod_U_t60_140.npz \
    --time-start 60 \
    --time-step 0.5 \
    --energy-threshold 0.9999 \
    --n-modes-max 30

python python/scripts/run_cnm_analysis.py \
    --input-type npz \
    --input data/processed/pod_U_t60_140.npz \
    --out data/figures/cnm_field_K10 \
    --K 10
```

For a sweep over time windows and cluster counts:

```bash
K_LIST="6 8 10 12" STRIDE=10 ./scripts/sweep_cnm_windows.sh
```

## Outputs

Each run writes a local result directory containing:

| File | Meaning |
|---|---|
| `01_pod_energy.png` | modal and cumulative POD energy |
| `02_mds_clusters.png` | MDS visualization |
| `03_P_matrix.png` | sampled CMM transition matrix |
| `04_Q_matrix.png` | direct CNM transition matrix |
| `05_T_matrix.png` | mean direct-transition times |
| `06_cluster_affiliation.png` | data and predicted cluster labels |
| `07_first_coefficients.png` | leading coefficients and reconstruction |
| `08_cluster_probabilities.png` | measured and predicted occupations |
| `summary.json` | run parameters and basic statistics |
| `*.npy` | numeric intermediate results |

These files are reproducible outputs and are ignored by Git.

## Interpretation limits

- The scripts currently train and reconstruct on the same time series. The
  plots are descriptive and are not an out-of-sample validation.
- Results depend on the chosen time window, stride, number of modes, and $K$.
- A hard POD mode cap can retain much less energy than the requested threshold.
- Euclidean k-means assumes that the selected coordinates and their scaling
  define a physically meaningful distance.
- The deterministic `argmax(Q)` trajectory suppresses less likely transitions.
- DMD and the importable `python/src/` library are placeholders, not finished
  implementations.

## Keeping Git small

The `.gitignore` is designed so that a normal `git add .` only stages source
and documentation. Before a commit, verify with:

```bash
git status --short
git check-ignore -v data/processed/example.npy
git check-ignore -v openfoam/cases/cylinderLES3D/VTK/example.vtk
```

A source-only ZIP can be created with:

```bash
./package_project.sh
```

The script writes the archive to `dist/`, which is itself ignored.

## References

Li, H., Fernex, D., Semaan, R., Tan, J., Morzynski, M., and Noack, B. R.
(2021). Cluster-based network model. *Journal of Fluid Mechanics*, 906, A21.
<https://doi.org/10.1017/jfm.2020.785>

Taira, K. et al. (2017). Modal analysis of fluid flows: An overview.
*AIAA Journal*, 55(12).
