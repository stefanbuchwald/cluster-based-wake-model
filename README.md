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
`-- postprocessing.sh
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

```math
\boldsymbol{x}(t_i)
=
\begin{bmatrix}
C_D(t_i) & C_L(t_i) & C_M(t_i)
\end{bmatrix}.
```

Because these observables have different scales, every column is standardized:

```math
\widetilde{x}_{ij}
=
\frac{x_{ij}-\mu_j}{s_j},
```

where $\mu_j$ and $s_j$ are the sample mean and standard deviation of feature
$j$.

### Velocity or pressure fields

`field_pod_from_vtk.py` reads cell data from a sequence of VTK snapshots.
Velocity components and optional pressure values are flattened into one
snapshot vector:

```math
\boldsymbol{x}_i
=
\begin{bmatrix}
U_x & U_y & U_z & p
\end{bmatrix}_{t_i}^{\mathsf T}.
```

The pressure block can be multiplied by `--p-scale`. This matters because a
combined velocity-pressure norm otherwise mixes quantities with different
units and physical meanings.

The complete snapshot matrix is stored locally as a memory-mapped `.npy` file.
It is intentionally ignored by Git because it can be several gigabytes large.

## Method

### POD and clustering

The input time series is centered and projected onto POD coordinates. Small
feature matrices are decomposed directly by SVD; large VTK fields use the
method of snapshots to avoid constructing a full spatial covariance matrix.
The retained temporal coefficients form the state space for clustering.

The number of retained modes is selected from the requested energy threshold,
subject to `--n-modes-max`. A hard mode limit can therefore stop before the
threshold is reached. For example, the existing 161-snapshot velocity-field
run retained about 89.27% of the fluctuation energy in its first 30 modes.

The reduced states are partitioned with k-means++ using 50 initializations and
a fixed random seed. Cluster numbers are then reordered for readability: the
most populated cluster becomes $C_0$, followed where possible by its most
frequent outgoing transition. This relabeling does not change the partition.

### Transition model

Let $l_i$ be the cluster label at sample $i$. For every pair of consecutive
samples, define

```math
n_{jk}
=
\#\left\{i:l_i=k,\ l_{i+1}=j\right\}.
```

The sampled transition matrix is

```math
P_{jk}
=
\frac{n_{jk}}{\sum_m n_{mk}}.
```

Thus, column $k$ contains the probability of the next sampled state, given
that the current state is $C_k$. Because consecutive samples in the same
cluster are counted, the diagonal of $P$ measures persistence at the chosen
sampling interval.

For the cluster-based network model, repeated labels are first compressed into
runs:

```text
0, 0, 0, 2, 2, 1, 1  ->  0, 2, 1
```

Let $\widehat n_{jk}$ count transitions between consecutive runs. Self
transitions no longer occur, and the direct-transition matrix becomes

```math
Q_{jk}
=
\frac{\widehat n_{jk}}
{\sum_{m\ne k}\widehat n_{mk}},
\qquad Q_{kk}=0.
```

Columns with no observed departure remain zero.

This distinction is important:

- `P` describes the next sampled label and depends directly on the sampling
  interval.
- `Q` describes the next different cluster after the current cluster is left.
- `T` supplies the characteristic time needed for that direct transition.

### Transition times

Let run $r$ start at sample $s_r$ and end at sample $e_r$. With

```math
\Delta t_{\mathrm{med}}
=
\mathrm{median}\left(t_{i+1}-t_i\right),
```

the implemented residence duration is

```math
R_r
=
t_{e_r}-t_{s_r}+\Delta t_{\mathrm{med}}.
```

For an observed run transition $k\rightarrow j$, the individual transition
time is

```math
\tau_r
=
\frac{R_r+R_{r+1}}{2}.
```

The matrix entry is the conditional average over all observed
$k\rightarrow j$ transitions:

```math
T_{jk}
=
\frac{1}{\widehat n_{jk}}
\sum_{r:\,k\rightarrow j}\tau_r.
```

Unobserved transitions are stored as `NaN`.

### Deterministic reconstruction

Starting from cluster $k$, the current implementation chooses the most likely
direct successor:

```math
j^\star
=
\arg\max_{j} Q_{jk}.
```

During the interval

```math
t_{\mathrm{arrive}}-t_{\mathrm{depart}}
=
T_{j^\star k},
```

the state is linearly interpolated between centroids:

```math
\widehat{\boldsymbol{a}}(t)
=
(1-\alpha)\boldsymbol{c}_k
+\alpha\boldsymbol{c}_{j^\star},
\qquad
\alpha
=
\frac{t-t_{\mathrm{depart}}}{T_{j^\star k}}.
```

The displayed cluster label switches from $k$ to $j^\star$ at
$\alpha=0.5$. If a cluster has no observed successor, the predicted state
remains at its centroid.

This reconstruction is deterministic. It does not sample from $Q$ and
therefore cannot reproduce stochastic branching or provide uncertainty
estimates.

### Visualization

Multidimensional scaling embeds snapshots and centroids into three dimensions
for the cluster map. It is only a visualization step and does not affect the
clustering, transition matrices, or reconstruction.

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

## Limitations and scope

This repository is an exploratory ROM workflow, not yet a predictive model
with quantified uncertainty.

- **No out-of-sample validation.** Clusters, transition statistics, and the
  displayed reconstruction are currently obtained from the same time series.
  The plots show how the fitted network represents its training data; they do
  not establish predictive performance on an unseen interval or another CFD
  run.
- **Dependence on temporal sampling.** The time window and `--stride` determine
  which transitions are observed. Coarse sampling can skip short residence
  events, while a window containing startup transients can produce a different
  network from a statistically developed wake.
- **Finite transition statistics.** Entries of `Q` and `T` are empirical
  estimates. Rare transitions may be based on very few observations, and a
  cluster with no observed departure becomes an absorbing state in the
  deterministic reconstruction. The current output does not report transition
  counts or confidence intervals.
- **State representation matters.** A model built from $C_D$, $C_L$, and $C_M$
  describes force dynamics, but distinct flow fields can share similar force
  coefficients. Field-based POD contains more spatial information, while its
  result depends on the retained modes and on any scaling between velocity and
  pressure.
- **Deterministic centroid dynamics.** The reconstruction always follows
  `argmax(Q)` and interpolates linearly between cluster centroids. It therefore
  suppresses less likely branches and does not reproduce within-cluster
  fluctuations, turbulent variability, or uncertainty.

Accordingly, the current results should be interpreted as a compact,
data-driven description of recurrent wake states and their observed
transitions. A predictive assessment would require separate training and test
windows, transition-count diagnostics, and comparison against withheld CFD
data.

## Keeping Git small

The `.gitignore` is designed so that a normal `git add .` only stages source
and documentation. Before a commit, verify with:

```bash
git status --short
git check-ignore -v data/processed/example.npy
git check-ignore -v openfoam/cases/cylinderLES3D/VTK/example.vtk
```

## References

Li, H., Fernex, D., Semaan, R., Tan, J., Morzynski, M., and Noack, B. R.
(2021). Cluster-based network model. *Journal of Fluid Mechanics*, 906, A21.
<https://doi.org/10.1017/jfm.2020.785>

Taira, K. et al. (2017). Modal analysis of fluid flows: An overview.
*AIAA Journal*, 55(12).
