#!/usr/bin/env pvpython
"""
export_lambda2_pvd.py

Robuster Export einer Lambda2-Zeitserie aus einem OpenFOAM-Case nach VTU + PVD.

Usage:
  pvpython export_lambda2_pvd.py case.foam [outdir] [stride] [--binary] [--stats]

Beispiele:
  pvpython export_lambda2_pvd.py case.foam lambda2_series 1 --binary
  pvpython export_lambda2_pvd.py case.foam lambda2_series 5 --stats
"""
from paraview.simple import *
import os, sys, math

args = [a for a in sys.argv[1:] if a != "--"]
case   = args[0] if len(args) > 0 else "case.foam"
outdir = args[1] if len(args) > 1 else "lambda2_series"
stride = int(args[2]) if len(args) > 2 and args[2].lstrip("-").isdigit() else 1

binary = ("--binary" in sys.argv[1:]) or ("-b" in sys.argv[1:])
stats  = ("--stats"  in sys.argv[1:]) or ("-s" in sys.argv[1:])

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

# --- Reader ---
r = OpenFOAMReader(FileName=case)
r.MeshRegions = ['internalMesh']

# U kann je nach Reader als Cell oder Point ankommen -> beides versuchen zu aktivieren
try:
    r.CellArrays = ['U']
except Exception:
    pass
try:
    r.PointArrays = ['U']
except Exception:
    pass

UpdatePipeline()

scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()

# --- Merge blocks (Multiblock -> UnstructuredGrid) ---
mb = MergeBlocks(Input=r)
UpdatePipeline()

# --- Gradient (CELLS -> fallback POINTS) ---
grad = Gradient(Input=mb)
grad.ResultArrayName = 'GradU'

ok = False
for assoc in ('CELLS', 'POINTS'):
    try:
        grad.ScalarArray = [assoc, 'U']
        UpdatePipeline()
        ok = True
        print("[export] Gradient uses", assoc, "U")
        break
    except Exception as e:
        print("[export] Gradient failed for", assoc, "U ->", e)

if not ok:
    raise RuntimeError("Could not compute gradient: U not found as CELLS or POINTS on merged dataset.")

# --- Lambda2 ProgrammableFilter (arbeitet mit CellData ODER PointData GradU) ---
pf = ProgrammableFilter(Input=grad)
pf.OutputDataSetType = 'Same as Input'
pf.Script = r"""
import numpy as np
from vtk.numpy_interface import dataset_adapter as dsa

inp = dsa.WrapDataObject(inputs[0].VTKObject)

# GradU kann in CellData oder PointData liegen (je nach Gradient-Input)
G = None
assoc = None
if 'GradU' in inp.CellData.keys():
    G = np.asarray(inp.CellData['GradU'])
    assoc = 'CELL'
elif 'GradU' in inp.PointData.keys():
    G = np.asarray(inp.PointData['GradU'])
    assoc = 'POINT'
else:
    raise RuntimeError("GradU not found on input (neither CellData nor PointData)")

G = np.squeeze(G)
if G.size == 0:
    raise RuntimeError("GradU is empty")

# GradU kommt häufig als (N,9) oder (N,3,3)
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

# Schreibe Lambda2 in dieselbe Datenassoziation wie GradU
if assoc == 'CELL':
    output.CellData.append(lambda2, 'Lambda2')
else:
    output.PointData.append(lambda2, 'Lambda2')
"""
UpdatePipeline()

# --- Timesteps bestimmen ---
tsteps = list(getattr(r, "TimestepValues", []))
if not tsteps:
    tsteps = [scene.AnimationTime]

# In ParaView-Skripten kann `max` überschrieben sein -> kein max() verwenden
step = stride if stride and stride > 1 else 1
tsteps = tsteps[::step]


print("[export] Found timesteps:", len(tsteps), "stride:", stride)
print("[export] Output:", os.path.abspath(outdir))

# --- optional stats file ---
stats_path = os.path.join(outdir, "lambda2_ranges.csv")
if stats:
    with open(stats_path, "w") as sf:
        sf.write("index,time,assoc,min,max\n")

# --- PVD schreiben ---
pvd_path = os.path.join(outdir, "lambda2_series.pvd")
mode = "Binary" if binary else "Ascii"

with open(pvd_path, "w") as f:
    f.write('<?xml version="1.0"?>\n')
    f.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n')
    f.write('  <Collection>\n')

    for i, t in enumerate(tsteps):
        # explizit auf Zeit updaten (verhindert "kleben" am ersten timestep)
        scene.AnimationTime = t
        UpdatePipeline(time=t, proxy=pf)

        fname = f"lambda2_{i:05d}.vtu"
        full  = os.path.join(outdir, fname)

        SaveData(full, proxy=pf, DataMode=mode)
        f.write(f'    <DataSet timestep="{t}" group="" part="0" file="{fname}"/>\n')

        if stats:
            di = pf.GetClientSideObject().GetOutputDataObject(0)
            cd = di.GetCellData().GetArray("Lambda2")
            pd = di.GetPointData().GetArray("Lambda2")
            if cd is not None:
                r0, r1 = cd.GetRange()
                assoc = "CELL"
            elif pd is not None:
                r0, r1 = pd.GetRange()
                assoc = "POINT"
            else:
                r0, r1, assoc = math.nan, math.nan, "NONE"
            with open(stats_path, "a") as sf:
                sf.write(f"{i},{t},{assoc},{r0},{r1}\n")

        if i % 10 == 0:
            print("[export] Saved", fname, "time", t)

    f.write('  </Collection>\n')
    f.write('</VTKFile>\n')

print("[export] DONE. Open in ParaView:", pvd_path)
if stats:
    print("[export] Stats:", stats_path)
