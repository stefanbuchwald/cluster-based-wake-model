# Synthetic quickstart

Run this example from the repository root after installing the Python
requirements. It needs no OpenFOAM installation or simulation files:

```bash
python scripts/generate_demo.py
python python/scripts/run_cnm_analysis.py \
    --input-type npz \
    --input data/processed/demo_periodic.npz \
    --out data/figures/demo \
    --K 8 \
    --n-runs 20
python python/scripts/run_dmd_analysis.py \
    --input data/processed/demo_periodic.npz \
    --out data/figures/demo
python python/scripts/plot_release_diagnostics.py --out data/figures/demo
```

The generator creates 241 uniformly spaced samples over `t=0...60` with
`dt=0.25`. Its four features are sine/cosine pairs at frequencies `0.2` and
`0.4`, with harmonic amplitude `0.3`. These are **synthetic signals**, not
simulation results or physical velocity fields. With `D=U_inf=1`, exact DMD
should identify `St=0.2` and `St=0.4` with near-zero growth rates and a relative
reconstruction error close to floating-point precision. The CNM run creates
clusters, transition matrices, cycle diagnostics, and reconstruction plots.
Cluster labels can differ between dependency versions.

All generated arrays, figures, and summaries stay under the ignored `data/`
directory. Use a different `--out` on the generator and corresponding input
arguments if you want to keep several examples.

The [reference results](../docs/results/README.md) come from the actual
OpenFOAM run and are separate from this demonstration.
