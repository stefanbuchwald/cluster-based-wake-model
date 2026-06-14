#!/usr/bin/env pvpython
"""
lambda2_render.py

Render eines einzelnen Lambda2-Frames.
- Input kann .foam oder .pvd/.vtu sein.
- Wenn Lambda2 nicht existiert, wird es aus U berechnet.

Usage:
  pvpython lambda2_render.py <case.foam|lambda2_series.pvd> [iso] [out.png] [--time <timestep>]

Beispiele:
  pvpython lambda2_render.py case.foam -1.0 lambda2.png
  pvpython lambda2_render.py lambda2_series/lambda2_series.pvd -0.5 lambda2.png --time 0.25
"""
from paraview.simple import *
import sys, os

args = sys.argv[1:]

inp = args[0] if len(args) > 0 else "case.foam"
iso = float(args[1]) if len(args) > 1 else -1.0
outpng = args[2] if len(args) > 2 else "lambda2.png"

time_val = None
if "--time" in args:
    i = args.index("--time")
    if i+1 < len(args):
        time_val = float(args[i+1])

def has_array(proxy, assoc, name):
    try:
        di = proxy.GetDataInformation()
        if assoc == "POINTS":
            ai = di.GetPointDataInformation().GetArrayInformation(name)
        else:
            ai = di.GetCellDataInformation().GetArrayInformation(name)
        return ai is not None
    except Exception:
        return False

if inp.lower().endswith(".foam"):
    r = OpenFOAMReader(FileName=inp)
    r.MeshRegions = ['internalMesh']
    try:
        r.CellArrays = ['U', 'Lambda2']
    except Exception:
        pass
    try:
        r.PointArrays = ['U', 'Lambda2']
    except Exception:
        pass
    src = r
else:
    src = OpenDataFile(inp)

UpdatePipeline()

scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()

tsteps = list(getattr(src, "TimestepValues", []))
if not tsteps:
    try:
        tsteps = list(scene.TimeKeeper.TimestepValues)
    except Exception:
        tsteps = [scene.AnimationTime]

if time_val is None:
    time_val = tsteps[0]

scene.AnimationTime = time_val

need_compute = not (has_array(src, "CELLS", "Lambda2") or has_array(src, "POINTS", "Lambda2"))
if need_compute:
    mb = MergeBlocks(Input=src)
    UpdatePipeline()

    grad = Gradient(Input=mb)
    grad.ResultArrayName = 'GradU'

    ok = False
    for assoc in ('CELLS', 'POINTS'):
        try:
            grad.ScalarArray = [assoc, 'U']
            UpdatePipeline()
            ok = True
            break
        except Exception:
            pass
    if not ok:
        raise RuntimeError("Could not compute gradient: U not found as CELLS or POINTS.")

    pf = ProgrammableFilter(Input=grad)
    pf.OutputDataSetType = 'Same as Input'
    pf.Script = r"""
import numpy as np
from vtk.numpy_interface import dataset_adapter as dsa

inp = dsa.WrapDataObject(inputs[0].VTKObject)

G = None
assoc = None
if 'GradU' in inp.CellData.keys():
    G = np.asarray(inp.CellData['GradU']); assoc='CELL'
elif 'GradU' in inp.PointData.keys():
    G = np.asarray(inp.PointData['GradU']); assoc='POINT'
else:
    raise RuntimeError("GradU not found on input (neither CellData nor PointData)")

G = np.squeeze(G)
if G.size == 0:
    raise RuntimeError("GradU is empty")

if G.ndim == 2 and G.shape[1] == 9:
    G = G.reshape((-1,3,3))
elif G.ndim == 3 and G.shape[1:] == (3,3):
    pass
else:
    raise RuntimeError(f"Unexpected GradU shape: {G.shape}")

S = 0.5*(G + np.transpose(G,(0,2,1)))
O = 0.5*(G - np.transpose(G,(0,2,1)))
A = np.einsum('nij,njk->nik', S, S) + np.einsum('nij,njk->nik', O, O)

w = np.linalg.eigvals(A).real
w.sort(axis=1)
lambda2 = w[:,1]

if assoc == 'CELL':
    output.CellData.append(lambda2, 'Lambda2')
else:
    output.PointData.append(lambda2, 'Lambda2')
"""
    UpdatePipeline(time=time_val, proxy=pf)
    data_for_contour = pf
else:
    data_for_contour = src

UpdatePipeline()

# Lambda2 als POINTS
if has_array(data_for_contour, "POINTS", "Lambda2"):
    l2_points = data_for_contour
else:
    l2_points = CellDatatoPointData(Input=data_for_contour)
    UpdatePipeline()

cont = Contour(Input=l2_points)
cont.ContourBy = ['POINTS', 'Lambda2']
cont.Isosurfaces = [iso]
UpdatePipeline(time=time_val, proxy=cont)

view = CreateView('RenderView')
view.ViewSize = [1600, 900]
disp = Show(cont, view)
disp.Representation = 'Surface'
ColorBy(disp, None)
ResetCamera(view)
view.ViewTime = time_val
Render()
SaveScreenshot(outpng, view)

print("Wrote", outpng, "| iso:", iso, "| time:", time_val)
