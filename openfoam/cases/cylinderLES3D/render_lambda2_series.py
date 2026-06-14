#!/usr/bin/env pvpython
"""
render_lambda2_series.py

Rendert eine Lambda2-Isofläche über alle Timesteps.
- Input kann .foam (OpenFOAM) ODER bereits exportierte .pvd/.vtu-Serie sein.
- Wenn Lambda2 im Input schon existiert, wird nicht erneut berechnet.

Usage:
  pvpython render_lambda2_series.py <case.foam|lambda2_series.pvd> [iso|-0.5,-1,-2] [outdir] [--clip] [--no-camera]

Beispiele:
  pvpython render_lambda2_series.py case.foam -1.0 frames_lambda2
  pvpython render_lambda2_series.py lambda2_series/lambda2_series.pvd -0.5 frames_lambda2
  pvpython render_lambda2_series.py case.foam "-0.5,-1,-2" frames_lambda2
"""
from paraview.simple import *
import os, sys

args = [a for a in sys.argv[1:] if a != "--"]
inp    = args[0] if len(args) > 0 else "case.foam"
iso_in = args[1] if len(args) > 1 else "-1.0"
outdir = args[2] if len(args) > 2 else "frames_lambda2"

do_clip     = ("--clip" in sys.argv[1:])
use_camera  = ("--no-camera" not in sys.argv[1:])

# iso kann Zahl oder CSV-Liste sein
isos = []
for part in str(iso_in).split(","):
    part = part.strip()
    if part:
        isos.append(float(part))
if not isos:
    isos = [-1.0]

os.makedirs(outdir, exist_ok=True)

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

# --- Source ---
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
    # .pvd/.vtu/... laden
    src = OpenDataFile(inp)

UpdatePipeline()

scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()

# --- Lambda2 vorhanden? ---
need_compute = not (has_array(src, "CELLS", "Lambda2") or has_array(src, "POINTS", "Lambda2"))

if need_compute:
    print("[render] Lambda2 not found -> computing from U")

    # Merge blocks (falls Multiblock)
    mb = MergeBlocks(Input=src)
    UpdatePipeline()

    # Gradient (CELLS -> POINTS fallback)
    grad = Gradient(Input=mb)
    grad.ResultArrayName = 'GradU'

    ok = False
    for assoc in ('CELLS', 'POINTS'):
        try:
            grad.ScalarArray = [assoc, 'U']
            UpdatePipeline()
            ok = True
            print("[render] Gradient uses", assoc, "U")
            break
        except Exception as e:
            print("[render] Gradient failed for", assoc, "U ->", e)
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
    data_for_contour = pf
else:
    print("[render] Lambda2 found on input -> using it directly")
    data_for_contour = src

# --- sicherstellen: Lambda2 als POINTS für Contour ---
UpdatePipeline()
if has_array(data_for_contour, "POINTS", "Lambda2"):
    l2_points = data_for_contour
else:
    # Cell->Point
    l2_points = CellDatatoPointData(Input=data_for_contour)
    UpdatePipeline()

# --- optional Clip ---
if do_clip:
    clip = Clip(Input=l2_points)
    clip.ClipType = 'Plane'
    clip.ClipType.Normal = [1, 0, 0]
    clip.ClipType.Origin = [0, 0, 1.0]
    UpdatePipeline()
    geom = clip
else:
    geom = l2_points

# --- Contour (mehrere Isos möglich) ---
cont = Contour(Input=geom)
cont.ContourBy = ['POINTS', 'Lambda2']
cont.Isosurfaces = isos
UpdatePipeline()

# --- Render setup ---
view = CreateView('RenderView')
view.ViewSize = [1920, 1080]
disp = Show(cont, view)
disp.Representation = 'Surface'
ColorBy(disp, None)

# Kamera (optional) – falls du aus GUI extrahiert hast
if use_camera:
    view.CameraParallelProjection = 0
    view.CameraPosition = [35.98050427826338, 22.41173854872093, 18.963597150478044]
    view.CameraFocalPoint = [13.028557619109293, -0.5402081104331762, -3.9883495086760714]
    view.CameraViewUp = [-0.4082482904638631, 0.816496580927726, -0.40824829046386296]
    view.CameraViewAngle = 35
else:
    ResetCamera(view)

Render()

# --- Timesteps ---
tsteps = list(getattr(src, "TimestepValues", []))
if not tsteps:
    # ParaView speichert Timesteps oft im TimeKeeper
    try:
        tsteps = list(scene.TimeKeeper.TimestepValues)
    except Exception:
        tsteps = [scene.AnimationTime]

print("[render] Timesteps:", len(tsteps), "| isos:", isos, "| out:", os.path.abspath(outdir))

for i, t in enumerate(tsteps):
    scene.AnimationTime = t
    view.ViewTime = t
    # explizit pipeline updaten – sonst bleibt manchmal frame 0
    UpdatePipeline(time=t, proxy=cont)
    Render()
    SaveScreenshot(os.path.join(outdir, f"frame_{i:05d}.png"), view)
    if i % 10 == 0:
        print("[render] Saved", i, "time", t)

print("[render] Done.")
