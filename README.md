# CFD wake reduced-order model

A small research and teaching project for building a cluster-based
reduced-order model (ROM) of a cylinder wake.

The repository contains source code and configuration only. Simulation fields,
VTK files, POD matrices, plots, and other generated results are deliberately
excluded from Git.

The workflow is:

```text
OpenFOAM LES, t = 0 ... 300
                  |
                  v
force coefficients and field snapshots
                  |
                  v
analysis window, t = 140 ... 270
                  |
                  v
        POD/SVD of velocity snapshots
                  |
                  v
           k-means++ clusters
                  |
                  v
     transition matrices P, Q, and T
                  |
                  v
 CNM cycles, probabilistic CMM, and DMD comparison
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
|       |-- export_lambda2_pvd.py
|       `-- lambda2_stats.py
|-- python/
|   |-- scripts/
|   |   |-- dmd_mode_field.py
|   |   |-- field_pod_from_vtk.py
|   |   |-- plot_force_coeffs.py
|   |   |-- run_dmd_analysis.py
|   |   `-- run_cnm_analysis.py
|   |-- src/rom/dmd.py
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
threshold is reached. For OpenFOAM field snapshots, `field_pod_from_vtk.py` uses
a finite-volume, volume-weighted POD by default. Schema-v2 field-POD NPZ files
already contain POD coefficients, so `run_cnm_analysis.py` truncates these
coefficients directly instead of applying a second POD. The current release run
uses `data/processed/pod_U_volume_t140_270.npz`, which keeps 200 modes and
retains about 99.71% of the volume-weighted fluctuation energy.

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

The deterministic reconstruction is retained as a reference. In addition, a
probabilistic CMM ensemble samples the full column $P_{\cdot k}$, including
self-transitions. Repeated realizations therefore retain stochastic branching
and provide a mean trajectory, a standard-deviation band, and an ensemble
cluster occupation. This is still an empirical in-sample model, not a calibrated
uncertainty estimate.

### Visualization

Multidimensional scaling embeds snapshots and centroids into three dimensions
for the cluster map. It is only a visualization step and does not affect the
clustering, transition matrices, or reconstruction.

### DMD comparison

The optional DMD runner applies exact Dynamic Mode Decomposition to the matrix
stored as `X` in the NPZ file. For files produced by
`field_pod_from_vtk.py`, this matrix contains retained temporal POD
coefficients, not the full spatial snapshot matrix. The resulting DMD modes
first describe dynamics in POD-coordinate space. Schema-v2 field POD files also
store the spatial POD basis and cell-center coordinates, allowing a selected
DMD mode to be projected back to the physical velocity or pressure field.

Discrete eigenvalues are converted to continuous growth rates and Strouhal
numbers. The runner reports independently ranked positive-frequency DMD
candidates; the release diagnostic plot then compares these candidates with
cycle frequencies measured from the supported CNM `Q/T` network.

## Running the analysis

From force coefficients:

```bash
python python/scripts/run_cnm_analysis.py \
    --input-type force \
    --input openfoam/cases/cylinderLES3D/postProcessing/forceCoeffs/0/forceCoeffs.dat \
    --out data/figures/cnm_force_K10_t140_270 \
    --K 10 \
    --t-min 140 \
    --t-max 270
```

From field POD coefficients, first export or keep VTK snapshots only for the
analysis window `t=140...270`; otherwise the glob will include the full local
VTK series:

```bash
python python/scripts/field_pod_from_vtk.py \
    --vtk-glob "openfoam/cases/cylinderLES3D/VTK/cylinderLES3D_*.vtk" \
    --fields U \
    --out data/processed/pod_U_volume_t140_270.npz \
    --time-start 140 \
    --time-step 0.25 \
    --energy-threshold 0.9999 \
    --n-modes-max 200 \
    --pod-weighting volume

python python/scripts/run_cnm_analysis.py \
    --input-type npz \
    --input data/processed/pod_U_volume_t140_270.npz \
    --out data/figures/cnm_field_K10 \
    --K 10 \
    --n-runs 200
```

For a sweep over time windows and cluster counts:

```bash
K_LIST="6 8 10 12" STRIDE=10 ./scripts/sweep_cnm_windows.sh
```

For an exact-DMD comparison using the field POD coefficients:

```bash
python python/scripts/run_dmd_analysis.py \
    --input data/processed/pod_U_volume_t140_270.npz \
    --out data/figures/dmd_U_volume_t140_270 \
    --dt 0.25

python python/scripts/dmd_mode_field.py \
    --npz data/processed/pod_U_volume_t140_270.npz \
    --dmd-summary data/figures/dmd_U_volume_t140_270/dmd_summary.json \
    --dmd-dir data/figures/dmd_U_volume_t140_270 \
    --mode-st 0.205403 \
    --component x \
    --z-slice 2.0 \
    --out data/figures/dmd_U_volume_t140_270/14_dmd_mode_field_St0205.png
```

The spatial plot requires a schema-v2 NPZ produced by the current
`field_pod_from_vtk.py`. Older NPZ files remain usable for CNM and
coefficient-space DMD, but must be regenerated once for physical DMD mode
plots.

For an interactive three-dimensional view, open the OpenFOAM case directly in
ParaView:

```bash
paraview openfoam/cases/cylinderLES3D/cylinderLES3D.foam
```

Select `Lambda2` and `U` in the OpenFOAM reader and create a contour at
$\lambda_2=-0.25$. Camera, lighting, coloring, and screenshot composition can
then be adjusted interactively. Data generation and the exact GUI steps are
documented in [`openfoam/README.md`](openfoam/README.md).

### Reproducing the final analysis

The final result discussed in the project overview uses the completed
OpenFOAM run from `t=0` to `t=300` and applies the ROM analysis to the
developed interval `t=140` to `t=270`:

```bash
source .venv/bin/activate
export MPLCONFIGDIR=/tmp/matplotlib

OUT="data/figures/final_release_volume_t140_270_K10"
NPZ="data/processed/pod_U_volume_t140_270.npz"
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
    --random-state 1 \
    --min-cycle-count 2

python python/scripts/run_dmd_analysis.py \
    --input "$NPZ" \
    --out "$OUT" \
    --dt 0.25

python python/scripts/plot_release_diagnostics.py \
    --out "$OUT"

# Optional spatial DMD slice for the selected mode:
python python/scripts/dmd_mode_field.py \
    --npz "$NPZ" \
    --dmd-summary "$OUT/dmd_summary.json" \
    --dmd-dir "$OUT" \
    --mode-st 0.205403 \
    --component x \
    --z-slice 2.0 \
    --out "$OUT/14_dmd_mode_field_St0205.png"
```

This reproduces the current reference folder without requiring fresh VTK
exports, as long as `pod_U_volume_t140_270.npz` is still present locally. The NPZ file
is a generated artifact and is intentionally not committed.

## Outputs

Each run writes a local result directory containing:

| File | Meaning |
|---|---|
| `01_pod_energy.png` | modal and cumulative POD energy |
| `02_mds_clusters.png` | MDS visualization |
| `03_P_matrix.png` | sampled CMM transition matrix |
| `04_Q_matrix.png` | direct CNM transition matrix |
| `05_T_matrix.png` | mean direct-transition times |
| `05b_transition_counts.png` | finite-sample support of direct transitions |
| `06_cluster_affiliation.png` | data and predicted cluster labels |
| `07_first_coefficients.png` | leading coefficients and reconstruction |
| `08_cluster_probabilities.png` | measured and predicted occupations |
| `09_dmd_eigenvalues.png` | discrete DMD eigenvalues and unit circle |
| `10_dmd_spectrum.png` | DMD amplitudes over Strouhal number |
| `11_dmd_partial_modal_reconstruction.png` | partial modal DMD reconstruction |
| `12_cluster_probabilities_probabilistic.png` | data, deterministic CNM, and probabilistic CMM occupations |
| `13_coefficients_probabilistic.png` | deterministic reconstruction and probabilistic ensemble band |
| `14_dmd_mode_field_St0205.png` | optional selected DMD mode projected to a physical field slice |
| `15_network_cycle_support.png` | dominant CNM cycle with transition support and times |
| `16_cnm_dmd_frequency_comparison.png` | compact comparison of independently obtained CNM and DMD frequencies |
| `cnm_summary.json` | CNM/CMM parameters, transition support, cycles, and provenance |
| `dmd_summary.json` | DMD rank, reconstruction diagnostics, candidate modes, and provenance |
| `network_cycles.json` | directed cycles found from the supported Q/T network |
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
  deterministic reconstruction. The corrected output reports direct transition
  counts so strong-looking probabilities can be interpreted with their sample
  support.
- **State representation matters.** A model built from $C_D$, $C_L$, and $C_M$
  describes force dynamics, but distinct flow fields can share similar force
  coefficients. Field-based POD contains more spatial information, while its
  result depends on the retained modes and on any scaling between velocity and
  pressure.
- **Centroid dynamics.** The deterministic reference follows `argmax(Q)`.
  The probabilistic CMM samples branches from $P$ with one fixed data time step
  per Markov step. It does not use CNM transition times $T$; persistence is
  already encoded by self-transitions in $P$. The ensemble spread represents
  empirical transition sampling, not complete physical or statistical
  uncertainty.
- **Single-operator DMD.** Exact DMD fits one autonomous linear operator over
  the complete analysis window. A trajectory whose oscillation envelope first
  grows and later decays is therefore not expected to be reconstructed
  accurately by one set of exponential modes. DMD amplitudes are fitted to the
  first snapshot rather than optimized over the complete time interval.

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

Tu, J. H. et al. (2014). On dynamic mode decomposition: Theory and
applications. *Journal of Computational Dynamics*, 1(2), 391-421.
