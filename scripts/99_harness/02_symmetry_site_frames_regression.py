#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path
from typing import Any


def _write_text(path: Path, s: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(s, encoding="utf-8")


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _run(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True)


def _load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _toy_cyclopropane_pdb() -> str:
    # Minimal PDB with 3 carbons in a triangle. RDKit infers bonds from CONECT.
    # Atom names intentionally match a fake "ph" cg mapping.
    return (
        "HETATM    1  C0  LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    2  C1  LIG A   1       1.500   0.000   0.000  1.00  0.00           C\n"
        "HETATM    3  C2  LIG A   1       0.500   1.300   0.000  1.00  0.00           C\n"
        "CONECT    1    2    3\n"
        "CONECT    2    1    3\n"
        "CONECT    3    1    2\n"
        "END\n"
    )


def main() -> None:
    ap = argparse.ArgumentParser(description="Regression: symmetric aromatic CGs must emit swapped site frames (CD1/CD2).")
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    ap.add_argument("--tag", type=str, default=str(int(time.time())))
    args = ap.parse_args()

    root = args.repo_root.resolve()
    outdir = root / "processed/99_harness" / f"symmetry_site_frames_{args.tag}"
    outdir.mkdir(parents=True, exist_ok=True)

    ligand_pdb = outdir / "LIG.pdb"
    cg_atommap = outdir / "cg_atommap.json"
    polar_sites = outdir / "polar_sites.json"
    out_json = outdir / "site_frames.json"

    _write_text(ligand_pdb, _toy_cyclopropane_pdb())
    _write_json(
        cg_atommap,
        {
            "0_0": {
                "cg": "ph",
                "lgd_sel": ["C0", "C1", "C2"],
                "correspond_resname": "PHE",
                "correspond_names": ["CG", "CD1", "CD2"],
            }
        },
    )
    _write_json(polar_sites, {"sites": []})

    build_frames_py = root / "scripts/02_polar_sites/03_build_ligand_site_frames.py"
    _run(
        [
            "uv",
            "run",
            "-p",
            "3.11",
            "python",
            str(build_frames_py),
            str(ligand_pdb),
            "--cg-atommap",
            str(cg_atommap),
            "--polar-sites",
            str(polar_sites),
            "-o",
            str(out_json),
        ]
    )

    payload = _load_json(out_json)
    sites = payload.get("sites", [])
    if not isinstance(sites, list):
        raise TypeError(f"Unexpected site_frames.json structure: sites is {type(sites)}")

    # We expect both the base frame and an additional symmetry-swapped frame (CD1<->CD2).
    if len(sites) != 2:
        raise AssertionError(f"Expected 2 site frames (base + swap) for ph; got {len(sites)}. Output: {out_json}")

    site_ids = [str(s.get("site_id")) for s in sites]
    if "cgmap:0_0" not in site_ids:
        raise AssertionError(f"Missing base site_id 'cgmap:0_0' in {site_ids}")

    want_swap = "cgmap:0_0:swap_CD1_CD2"
    if want_swap not in site_ids:
        raise AssertionError(f"Missing swapped site_id '{want_swap}' in {site_ids}")

    # Ensure swapped frame uses the opposite CD1/CD2 assignment.
    swapped = next(s for s in sites if str(s.get("site_id")) == want_swap)
    frame_atoms = list(swapped.get("frame_atom_names", []))
    if frame_atoms != ["C0", "C2", "C1"]:
        raise AssertionError(f"Unexpected swapped frame_atom_names: {frame_atoms} (expected ['C0','C2','C1'])")


if __name__ == "__main__":
    main()

