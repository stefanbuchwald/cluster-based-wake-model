#!/usr/bin/env python3
"""Calculate (and optionally set) Reynolds number for an OpenFOAM case.

Designed for OpenFOAM.org v13-style dictionaries (physicalProperties) but also
works with older transportProperties.

Re = U_ref * D / nu

- U_ref: magnitude of inlet patch fixedValue (default patch: inlet) from 0/U
- D: cylinder diameter from a .geo file variable 'D = ...;' (best effort)
- nu: kinematic viscosity from constant/physicalProperties (or transportProperties)

Examples
--------
# just print current values
./calcRe.py

# show what you'd need for Re=200
./calcRe.py --target 200

# set nu so that Re becomes 200 (keep U and D)
./calcRe.py --target 200 --set-target-nu

# set inlet speed so that Re becomes 200 (keep nu and D)
./calcRe.py --target 200 --set-target-U

# directly set nu or inlet magnitude
./calcRe.py --set-nu 5e-3
./calcRe.py --set-U 0.05
"""

from __future__ import annotations

import argparse
import math
import os
import re
from pathlib import Path
from typing import Optional, Tuple


FLOAT_RE = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"


def read_text(p: Path) -> str:
    return p.read_text(encoding="utf-8", errors="ignore")


def find_case_root(case_dir: Path) -> Path:
    # Heuristic: must contain 'system' and 'constant'
    if (case_dir / "system").is_dir() and (case_dir / "constant").is_dir():
        return case_dir
    raise SystemExit(f"Not an OpenFOAM case directory: {case_dir}")


def parse_nu(case: Path) -> Tuple[float, Path]:
    """Return (nu, file_path)."""
    candidates = [
        case / "constant" / "physicalProperties",
        case / "constant" / "transportProperties",
    ]
    for p in candidates:
        if p.is_file():
            txt = read_text(p)
            m = re.search(rf"^\s*nu\s+({FLOAT_RE})\s*;", txt, flags=re.M)
            if m:
                return float(m.group(1)), p
            # Sometimes: nu [0 2 -1 0 0 0 0] 1e-6;
            m = re.search(rf"^\s*nu\s+\[[^\]]+\]\s+({FLOAT_RE})\s*;", txt, flags=re.M)
            if m:
                return float(m.group(1)), p
    raise SystemExit("Could not find 'nu' in constant/physicalProperties or constant/transportProperties")


def extract_patch_block(text: str, patch: str) -> Optional[str]:
    """Extract a patch sub-dictionary block like 'inlet { ... }' (best effort)."""
    lines = text.splitlines()
    in_block = False
    brace = 0
    buf = []

    # start must match: optional spaces, patch name, optional spaces, '{'
    start_pat = re.compile(rf"^\s*{re.escape(patch)}\s*\{{\s*$")

    for ln in lines:
        if not in_block:
            if start_pat.match(ln):
                in_block = True
                brace = ln.count("{") - ln.count("}")
                buf.append(ln)
        else:
            buf.append(ln)
            brace += ln.count("{") - ln.count("}")
            if brace == 0:
                return "\n".join(buf)

    return None


def parse_inlet_U(case: Path, patch: str = "inlet") -> Tuple[float, Tuple[float, float, float]]:
    """Return (|U|, (Ux,Uy,Uz)) from 0/U fixedValue uniform vector."""
    Ufile = case / "0" / "U"
    if not Ufile.is_file():
        raise SystemExit(f"Missing {Ufile}")

    txt = read_text(Ufile)

    block = extract_patch_block(txt, patch)
    if block is None:
        # fallback: search anywhere for 'patch { ... value uniform (..); }'
        block = txt

    m = re.search(
        rf"\bvalue\s+uniform\s*\(\s*({FLOAT_RE})\s+({FLOAT_RE})\s+({FLOAT_RE})\s*\)",
        block,
    )
    if not m:
        raise SystemExit(
            f"Could not find 'value uniform (Ux Uy Uz)' for patch '{patch}' in 0/U"
        )

    ux, uy, uz = (float(m.group(1)), float(m.group(2)), float(m.group(3)))
    mag = math.sqrt(ux * ux + uy * uy + uz * uz)
    return mag, (ux, uy, uz)


def parse_D_from_geo(case: Path) -> Optional[Tuple[float, Path]]:
    """Try to find cylinder diameter D from a *.geo file (Gmsh).

    Looks for a line like: D = 1.0;
    """
    geos = sorted(case.glob("*.geo"))
    for p in geos:
        txt = read_text(p)
        m = re.search(rf"^\s*D\s*=\s*({FLOAT_RE})\s*;", txt, flags=re.M)
        if m:
            return float(m.group(1)), p
    return None


def replace_first(pattern: str, repl: str, text: str) -> Tuple[str, bool]:
    new_text, n = re.subn(pattern, repl, text, count=1, flags=re.M)
    return new_text, (n == 1)


def set_nu_in_file(nu_file: Path, new_nu: float) -> None:
    txt = read_text(nu_file)
    # keep original formatting as much as possible
    new_txt, ok = replace_first(
        rf"^(\s*nu\s+)(?:\[[^\]]+\]\s+)?({FLOAT_RE})(\s*;)",
        rf"\g<1>{new_nu:.15g}\g<3>",
        txt,
    )
    if not ok:
        raise SystemExit(f"Could not update nu in {nu_file}")
    nu_file.write_text(new_txt, encoding="utf-8")


def set_inlet_U_magnitude(case: Path, patch: str, new_mag: float) -> None:
    Ufile = case / "0" / "U"
    txt = read_text(Ufile)

    # Extract current vector to preserve direction
    old_mag, (ux, uy, uz) = parse_inlet_U(case, patch=patch)

    if old_mag > 0:
        dx, dy, dz = ux / old_mag, uy / old_mag, uz / old_mag
    else:
        # default direction if zero
        dx, dy, dz = 1.0, 0.0, 0.0

    nx, ny, nz = new_mag * dx, new_mag * dy, new_mag * dz

    block = extract_patch_block(txt, patch)
    if block is None:
        # fallback: replace first occurrence in file
        new_txt, ok = replace_first(
            rf"(\bvalue\s+uniform\s*\()\s*{FLOAT_RE}\s+{FLOAT_RE}\s+{FLOAT_RE}(\s*\))",
            rf"\g<1>{nx:.15g} {ny:.15g} {nz:.15g}\g<2>",
            txt,
        )
        if not ok:
            raise SystemExit(f"Could not update inlet U in {Ufile}")
        Ufile.write_text(new_txt, encoding="utf-8")
        return

    # replace within that patch block only
    block_new, ok = replace_first(
        rf"(\bvalue\s+uniform\s*\()\s*{FLOAT_RE}\s+{FLOAT_RE}\s+{FLOAT_RE}(\s*\))",
        rf"\g<1>{nx:.15g} {ny:.15g} {nz:.15g}\g<2>",
        block,
    )
    if not ok:
        raise SystemExit(f"Could not update 'value uniform (...)' inside patch '{patch}'")

    # splice back
    txt_new = txt.replace(block, block_new)
    Ufile.write_text(txt_new, encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser(description="Compute Reynolds number for an OpenFOAM case")
    ap.add_argument("--case", default=".", help="case directory (default: .)")
    ap.add_argument("--patch", default="inlet", help="inlet patch name in 0/U (default: inlet)")
    ap.add_argument("--D", type=float, default=None, help="override cylinder diameter")
    ap.add_argument("--target", type=float, default=None, help="target Reynolds number")

    ap.add_argument("--set-nu", type=float, default=None, help="write nu to constant/physicalProperties")
    ap.add_argument("--set-U", type=float, default=None, help="write inlet |U| (keeps direction) to 0/U")

    ap.add_argument("--set-target-nu", action="store_true", help="set nu so Re hits --target (keeps U)")
    ap.add_argument("--set-target-U", action="store_true", help="set inlet |U| so Re hits --target (keeps nu)")

    args = ap.parse_args()

    case = find_case_root(Path(args.case).resolve())

    nu, nu_file = parse_nu(case)
    Uref, Uvec = parse_inlet_U(case, patch=args.patch)

    D = args.D
    D_src = None
    if D is None:
        got = parse_D_from_geo(case)
        if got:
            D, D_src = got

    if D is None:
        raise SystemExit(
            "Could not auto-detect D (cylinder diameter). Provide it via --D <value>."
        )

    # apply direct writes first
    if args.set_nu is not None:
        set_nu_in_file(nu_file, args.set_nu)
        nu = args.set_nu

    if args.set_U is not None:
        set_inlet_U_magnitude(case, args.patch, args.set_U)
        Uref = args.set_U

    # apply target-based writes
    if (args.set_target_nu or args.set_target_U) and args.target is None:
        raise SystemExit("--set-target-nu/--set-target-U requires --target")

    if args.set_target_nu and args.target is not None:
        new_nu = Uref * D / args.target
        set_nu_in_file(nu_file, new_nu)
        nu = new_nu

    if args.set_target_U and args.target is not None:
        new_U = args.target * nu / D
        set_inlet_U_magnitude(case, args.patch, new_U)
        Uref = new_U

    Re = Uref * D / nu

    print("=== Reynolds number (based on inlet |U|, cylinder D, molecular nu) ===")
    print(f"case     : {case}")
    print(f"patch    : {args.patch}")
    print(f"U_inlet  : {Uref:.15g}  (vector in file ~ {Uvec[0]:.6g} {Uvec[1]:.6g} {Uvec[2]:.6g})")
    print(f"nu       : {nu:.15g}  (from {nu_file.relative_to(case)})")
    if D_src:
        print(f"D        : {D:.15g}  (from {D_src.name})")
    else:
        print(f"D        : {D:.15g}")
    print(f"Re       : {Re:.6g}")

    if args.target is not None:
        nu_need = Uref * D / args.target
        U_need = args.target * nu / D
        print("\n--- Target helper ---")
        print(f"target Re: {args.target:.6g}")
        print(f"nu needed (keep U,D): {nu_need:.15g}")
        print(f"U needed  (keep nu,D): {U_need:.15g}")


if __name__ == "__main__":
    main()
