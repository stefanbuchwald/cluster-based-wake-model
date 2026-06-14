# OpenFOAM cylinder LES

This directory contains the CFD part of the project: a three-dimensional,
incompressible large-eddy simulation of flow past a circular cylinder.

Generated meshes, time directories, decomposed processor data, VTK exports,
and post-processing results are excluded from Git. The committed files are
enough to regenerate the case.

## Physical problem

The simulation uses the filtered incompressible Navier-Stokes equations:

$$
\nabla\cdot\overline{\boldsymbol{u}}=0,
$$

$$
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
$$

with strain-rate tensor

$$
\overline{\mathbf S}
=
\frac{1}{2}
\left(
\nabla\overline{\boldsymbol{u}}
+\nabla\overline{\boldsymbol{u}}^{\mathsf T}
\right).
$$

The reference quantities are

$$
U_\infty=1,\qquad D=1,\qquad
\nu=\frac{1}{3900}\approx2.5641\times10^{-4},
$$

and therefore

$$
\mathrm{Re}_D
=
\frac{U_\infty D}{\nu}
\approx3900.
$$

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

$$
\Delta=V_{\mathrm{cell}}^{1/3}.
$$

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
| `inlet` | fixed $(1,0,10^{-6})$ | zero gradient | uniform inflow with a tiny spanwise perturbation |
| `outlet` | zero gradient | fixed $0$ | pressure reference and convective outflow approximation |
| `walls` | slip | zero gradient | lateral far-field boundaries |
| `cylinder` | no slip | zero gradient | solid cylinder |
| `front`, `back` | cyclic | cyclic | periodic spanwise direction |

The spanwise perturbation is small enough not to alter the reference speed but
helps break an exactly two-dimensional initial symmetry.

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
| write interval | 0.5 |
| end time | 140 |

The initial `deltaT` is `0.0005`; afterwards OpenFOAM adapts it subject to the
Courant limit.

## Force coefficients

`system/controlDict` evaluates forces on the `cylinder` patch. With
$U_\infty=1$, $\rho=1$, $D=1$, and span $L_z=4$, the reference area is

$$
A_{\mathrm{ref}}=DL_z=4.
$$

The reported drag and lift coefficients follow

$$
C_D
=
\frac{F_x}{\tfrac12\rho U_\infty^2 A_{\mathrm{ref}}},
\qquad
C_L
=
\frac{F_y}{\tfrac12\rho U_\infty^2 A_{\mathrm{ref}}}.
$$

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
./calcRe.py
./calcRe.py --target 3900
```

Options that modify `nu` or the inlet velocity should be used deliberately,
because they change committed case configuration files.

## Post-processing

From the project root:

```bash
./postprocessing.sh --t-min 60 --K 10
```

This can reconstruct parallel results, export `U` and `p` to VTK, merge force
coefficient files, and run the Python CNM analysis.

Useful switches include:

```bash
./postprocessing.sh --no-reconstruct --no-vtk --t-min 60 --K 10
./postprocessing.sh --vtk-latest
./postprocessing.sh --vtk-time 60:140
```

For direct visualization:

```bash
cd openfoam/cases/cylinderLES3D
./Allview latest
```

`Allview` reconstructs the latest parallel time when processor directories
exist and then starts `paraFoam`.

The Lambda2 helper scripts in the case directory are optional ParaView/Python
tools. For vortex visualization, negative Lambda2 isosurfaces are commonly
used. Their generated `.vtu`, `.pvd`, frame, and video files are ignored.

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
    `-- optional Python post-processing helpers
```
