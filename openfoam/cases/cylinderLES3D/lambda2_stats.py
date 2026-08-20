#!/usr/bin/env pvpython
"""
lambda2_stats.py

Gibt (min,max) von Lambda2 aus.
- Input kann .foam oder .pvd/.vtu sein.
- Wenn Lambda2 nicht existiert, wird es aus U berechnet (Gradient + ProgrammableFilter).

Usage:
  pvpython lambda2_stats.py <case.foam|lambda2_series.pvd> [--per-time]

Beispiele:
  pvpython lambda2_stats.py case.foam
  pvpython lambda2_stats.py case.foam --per-time
  pvpython lambda2_stats.py lambda2_series/lambda2_series.pvd --per-time
"""
from paraview.simple import *
import sys, math

inp = sys.argv[1] if len(sys.argv) > 1 else "case.foam"
per_time = ("--per-time" in sys.argv[1:]) or ("-t" in sys.argv[1:])

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
    UpdatePipeline()
    l2 = pf
else:
    l2 = src

def range_at_time(t):
    scene.AnimationTime = t
    UpdatePipeline(time=t, proxy=l2)
    di = l2.GetClientSideObject().GetOutputDataObject(0)

    cd = di.GetCellData().GetArray("Lambda2")
    pd = di.GetPointData().GetArray("Lambda2")

    if cd is not None:
        r0, r1 = cd.GetRange()
        return "CELL", r0, r1
    if pd is not None:
        r0, r1 = pd.GetRange()
        return "POINT", r0, r1
    return "NONE", math.nan, math.nan

tsteps = list(getattr(src, "TimestepValues", []))
if not tsteps:
    try:
        tsteps = list(scene.TimeKeeper.TimestepValues)
    except Exception:
        tsteps = [scene.AnimationTime]

if not per_time:
    assoc, r0, r1 = range_at_time(tsteps[0])
    print("Lambda2 range:", (r0, r1), "| assoc:", assoc, "| time:", tsteps[0])
else:
    print("index,time,assoc,min,max")
    for i,t in enumerate(tsteps):
        assoc, r0, r1 = range_at_time(t)
        print(f"{i},{t},{assoc},{r0},{r1}")
