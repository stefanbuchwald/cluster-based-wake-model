# OpenFOAM cylinder LES

This directory contains the CFD part of the project: a three-dimensional,
incompressible large-eddy simulation of flow past a circular cylinder.

Generated meshes, time directories, decomposed processor data, VTK exports,
and post-processing results are excluded from Git. The committed files are
enough to regenerate the case.

## Physical problem

The simulation uses the filtered incompressible Navier-Stokes equations:

```math
\nabla\cdot\overline{\boldsymbol{u}}=0,
```

```math
\frac{\partial\overline{\boldsymbol{u}}}{\partial t}
+\nabla\cdot
\left(
\overline{\boldsymbol{u}}\otimes\overline{\boldsymbol{u}}
\right)
=
-\nabla\overline p
+\nabla\cdot
\left[
2(\nu+\nu_{\mathrm{sgs}})
\overline{\mathbf S}
\right],
```

with strain-rate tensor

```math
\overline{\mathbf S}
=
\frac{1}{2}
\left(
\nabla\overline{\boldsymbol{u}}
+\nabla\overline{\boldsymbol{u}}^{\mathsf T}
\right).
```

The reference quantities are

```math
U_\infty=1,\qquad D=1,\qquad
\nu=\frac{1}{3900}\approx2.5641\times10^{-4},
```

and therefore

```math
\mathrm{Re}_D
=
\frac{U_\infty D}{\nu}
\approx3900.
```

OpenFOAM uses kinematic pressure $p/\rho$. Density is set to one for the force
coefficients.

## LES model

`constant/momentumTransport` selects the WALE subgrid-scale model:

```text
simulationType  LES;
model           WALE;
delta           cubeRootVol;
Cw              0.325;
```

WALE constructs an eddy viscosity from the local velocity-gradient tensor and
has the appropriate near-wall scaling without requiring a dynamic averaging
procedure. The filter width is

```math
\Delta=V_{\mathrm{cell}}^{1/3}.
```

## Domain and mesh

The Gmsh geometry is
`cases/cylinderLES3D/karman_cylinder_3D_structured.geo`.

| Quantity | Value |
|---|---:|
| cylinder diameter | $D=1$ |
| inlet position | $x=-10D$ |
| outlet position | $x=25D$ |
| channel height | $10D$ |
| spanwise length | $L_z=4D$ |
| near-cylinder block half-width | $2.5D$ |

The two-dimensional multi-block grid consists of an O-grid-like region around
the cylinder and structured outer blocks. It is extruded in the spanwise
direction to form hexahedra. Transfinite curve progressions refine the mesh
near the cylinder and in the near wake.

`TargetCellsM` near the top of the `.geo` file controls the approximate mesh
size. The committed value is `0.25`; the comments in that file describe the
scaling.

`Allbuild` performs:

1. Gmsh mesh generation in MSH2 format.
2. Renumbering of physical tags to avoid surface/volume tag collisions.
3. Conversion with `gmshToFoam`.
4. Conversion of `front` and `back` to translational cyclic patches.
5. `checkMesh -allTopology -allGeometry`.

The generated `.msh` files and `constant/polyMesh/` are ignored by Git.

## Boundary conditions

| Patch | Velocity | Pressure | Meaning |
|---|---|---|---|
| `inlet` | fixed $(1,0,10^{-6})$ | zero gradient | uniform inflow with a tiny spanwise component |
| `outlet` | zero gradient | fixed $0$ | pressure reference and convective outflow approximation |
| `walls` | slip | zero gradient | lateral far-field boundaries |
| `cylinder` | no slip | zero gradient | solid cylinder |
| `front`, `back` | cyclic | cyclic | periodic spanwise direction |

The spanwise component is small enough not to alter the reference speed. It is
uniform in space and does not itself break spanwise translation symmetry. The
computed wake must be checked for spanwise variation when assessing the
three-dimensional flow.

## Numerical method

The active configuration is in `cases/cylinderLES3D/system/`.

| Setting | Value |
|---|---|
| OpenFOAM | OpenFOAM.org v13 |
| solver module | `incompressibleFluid` through `foamRun` |
| pressure-velocity coupling | PIMPLE |
| outer correctors | 2 |
| pressure correctors | 2 |
| time discretization | first-order Euler |
| convection of velocity | `linearUpwindV grad(U)` |
| gradients/interpolation | linear |
| Laplacian | linear corrected |
| adaptive time step | yes |
| maximum Courant number | 0.75 |
| maximum time step | 0.01 |
| write interval | 0.25 |
| end time | 300 |

The initial `deltaT` is `0.0005`; afterwards OpenFOAM adapts it subject to the
Courant limit.

## Force coefficients

`system/controlDict` evaluates forces on the `cylinder` patch. With
$U_\infty=1$, $\rho=1$, $D=1$, and span $L_z=4$, the reference area is

```math
A_{\mathrm{ref}}=DL_z=4.
```

The reported drag and lift coefficients follow

```math
C_D
=
\frac{F_x}{\tfrac12\rho U_\infty^2 A_{\mathrm{ref}}},
\qquad
C_L
=
\frac{F_y}{\tfrac12\rho U_\infty^2 A_{\mathrm{ref}}}.
```

The coefficient files are written below
`postProcessing/forceCoeffs/` and are ignored by Git.

## Requirements

- OpenFOAM.org v13
- Gmsh
- MPI for parallel runs
- Python 3 for the mesh-tag helper in `Allbuild`
- ParaView or `paraFoam` for visualization

Load OpenFOAM before running the case:

```bash
source /opt/openfoam13/etc/bashrc
```

Adjust the path if OpenFOAM is installed elsewhere.

## Build and run

```bash
cd openfoam/cases/cylinderLES3D

# Create the mesh and run checkMesh.
./Allbuild

# Serial run. Allrun rebuilds the mesh first.
./Allrun

# Parallel run with four MPI ranks.
./Allrun 4
```

For a parallel run, `Allrun` updates `numberOfSubdomains`, runs
`decomposePar -force`, executes `foamRun -parallel`, and reconstructs the
latest time after completion.

The Reynolds-number helper can inspect or modify the dimensional parameters:

```bash
python3 calcRe.py
python3 calcRe.py --target 3900
```

Options that modify `nu` or the inlet velocity should be used deliberately,
because they change committed case configuration files.

## Post-processing

From the project root:

```bash
./postprocessing.sh --t-min 140 --t-max 270 --K 10 --resample-dt 0.01 --tag final_t140_270
```

This can reconstruct parallel results, export `U` and `p` to VTK, merge force
coefficient files, and run the Python CNM analysis. Reconstruction is skipped
automatically for serial cases without processor directories.

Adaptive solver time steps also produce irregular force-coefficient times.
Choose `--resample-dt 0.01` to linearly interpolate them onto a fixed grid before
`--stride` is applied; the default stride of 10 then gives `dt = 0.1`. Omit this
option only for already uniform force data. Use `./postprocessing.sh --help` for
all options.

Useful switches include:

```bash
./postprocessing.sh --no-reconstruct --no-vtk --t-min 140 --t-max 270 --K 10 --resample-dt 0.01
./postprocessing.sh --vtk-latest --resample-dt 0.01
./postprocessing.sh --vtk-time 140:270 --resample-dt 0.01
```

### Lambda2 vortex structures

The Lambda2 criterion identifies vortex cores from the second eigenvalue of

```math
\mathbf S^2+\boldsymbol{\Omega}^2,
```

where $\mathbf S$ and $\boldsymbol{\Omega}$ are the symmetric and
antisymmetric parts of the velocity-gradient tensor. Connected regions with
$\lambda_2<0$ indicate rotation-dominated flow structures.

The case may contain local `Lambda2` fields in generated time directories.
They are ignored by Git. To regenerate them from the velocity field for the
developed interval:

```bash
cd openfoam/cases/cylinderLES3D

postProcess -func Lambda2 -time '140:270'
touch cylinderLES3D.foam
```

Open the case directly:

```bash
paraview cylinderLES3D.foam
```

In ParaView:

1. Select `internalMesh`, `U`, and `Lambda2` in the OpenFOAM reader and click
   **Apply**.
2. Move to a developed time, for example `t = 200`.
3. Add **Cell Data to Point Data** if `Lambda2` is available only as cell data.
4. Add a **Contour**, select `Lambda2`, and start with an isovalue of `-0.25`.
5. Color the surface by `U` using **Magnitude**, then adjust camera, lighting,
   background, and material interactively.
6. Save the final view with **File > Save Screenshot**.

For this case, `-0.25` retains the alternating wake structures over a
substantial downstream distance. More negative thresholds such as `-0.5` or
`-1.0` isolate progressively stronger and more compact vortex cores.

The available value range can be checked before opening ParaView:

```bash
pvpython lambda2_stats.py cylinderLES3D.foam
pvpython lambda2_stats.py cylinderLES3D.foam --per-time
```

ParaView can read the `.foam` time series directly. A VTU/PVD export is only
needed when the visualization data should be portable without the OpenFOAM
case:

```bash
pvpython export_lambda2_pvd.py \
    cylinderLES3D.foam \
    lambda2_series \
    5 \
    --binary \
    --stats
```

For direct visualization:

```bash
cd openfoam/cases/cylinderLES3D
./Allview latest
```

`Allview` reconstructs the latest parallel time when processor directories
exist and then starts `paraFoam`.

To inspect or stop solver processes belonging to this case:

```bash
./Allkill --dry-run
./Allkill
```

`Allkill` uses the explicit `-case` argument or the process working directory
and only selects OpenFOAM solvers and their MPI launchers owned by the current
user. It sends SIGTERM, then SIGKILL to matching processes still running after
one second. It requires Linux `/proc` and Python 3.

## Cleaning local results

`Allclean` removes the generated mesh, processor directories, numerical time
directories, MSH files, and logs:

```bash
./Allclean
```

This is destructive and should only be run when the local simulation results
are no longer needed. The Git ignore rules do not delete any data; they only
prevent generated files from being committed.

## Files kept in Git

```text
openfoam/
|-- README.md
`-- cases/cylinderLES3D/
    |-- 0/U
    |-- 0/p
    |-- 0/nut
    |-- constant/physicalProperties
    |-- constant/momentumTransport
    |-- system/controlDict
    |-- system/fvSchemes
    |-- system/fvSolution
    |-- system/decomposeParDict
    |-- system/createPatchDict
    |-- karman_cylinder_3D_structured.geo
    |-- Allbuild
    |-- Allrun
    |-- Allclean
    |-- Allkill
    |-- Allview
    |-- calcRe.py
    |-- export_lambda2_pvd.py
    `-- lambda2_stats.py
```
