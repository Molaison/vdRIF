#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path
from typing import Any

import numpy as np


def _run(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True)


def _load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    ap = argparse.ArgumentParser(description="Smoke: MTX candidates -> rif export inputs (placements + irot lib).")
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    ap.add_argument("--tag", type=str, default=str(int(time.time())))
    ap.add_argument("--top-per-site", type=int, default=200)
    ap.add_argument("--chunk-size", type=int, default=5000)
    args = ap.parse_args()

    root = args.repo_root.resolve()
    outdir = root / "processed/99_harness" / f"mtx_rif_export_{args.tag}"
    outdir.mkdir(parents=True, exist_ok=True)

    # Ensure inputs exist (polar + frames can be regenerated; vdxform must exist already).
    ligand = root / "inputs/01_cgmap/MTX.pdb"
    polar = root / "outputs/02_polar_sites/MTX_polar_sites.json"
    frames = root / "outputs/02_polar_sites/MTX_site_frames.json"
    vdx = root / "processed/03_vdxform"

    if not ligand.exists():
        raise FileNotFoundError(f"Missing {ligand}")
    for cg in ["coo", "conh2", "ccn", "ph", "bb_cnh"]:
        p = vdx / cg / f"vdxform_{cg}.npz"
        if not p.exists():
            raise FileNotFoundError(f"Missing {p}")

    # Regenerate polar-sites + frames (this script auto-falls back to rdkit when openbabel is unavailable).
    _run(["bash", str(root / "scripts/02_polar_sites/01_run_mtx_polar_sites.sh")])
    if not polar.exists() or not frames.exists():
        raise FileNotFoundError("Expected MTX polar-sites + site-frames to exist after 01_run_mtx_polar_sites.sh")

    cand_prefix = outdir / "MTX_candidates"
    candidate_py = root / "scripts/04_candidates/01_generate_candidates.py"
    export_py = root / "scripts/06_rif_export/01_export_rif_inputs.py"

    _run(
        [
            "uv",
            "run",
            "-p",
            "3.11",
            "python",
            str(candidate_py),
            "--ligand-pdb",
            str(ligand),
            "--polar-sites",
            str(polar),
            "--site-frames",
            str(frames),
            "--vdxform-dir",
            str(vdx),
            "--out-prefix",
            str(cand_prefix),
            "--chunk-size",
            str(int(args.chunk_size)),
            "--top-per-site",
            str(int(args.top_per_site)),
            "--clash-tol",
            "0.5",
            "--exclude-aa3",
            "PRO,CYS",
            "--require-sidechain-facing",
            "--min-sidechain-facing-dot",
            "0.2",
            "--require-full-coverage",
        ]
    )

    candidates_npz = cand_prefix.with_suffix(".npz")
    candidates_meta = cand_prefix.with_suffix(".json")
    if not candidates_npz.exists() or not candidates_meta.exists():
        raise FileNotFoundError("Candidate generation did not produce expected outputs.")

    out_prefix = outdir / "MTX"
    _run(
        [
            "uv",
            "run",
            "-p",
            "3.11",
            "python",
            str(export_py),
            "--candidates-npz",
            str(candidates_npz),
            "--candidates-meta",
            str(candidates_meta),
            "--site-frames",
            str(frames),
            "--out-prefix",
            str(out_prefix),
            "--rotamer-bits",
            "12",
            "--score-sign",
            "negate",
        ]
    )

    placements_npz = out_prefix.with_name(out_prefix.name + "_placements.npz")
    lib_npz = out_prefix.with_name(out_prefix.name + "_irot_lib.npz")
    meta_json = out_prefix.with_name(out_prefix.name + "_meta.json")
    for p in [placements_npz, lib_npz, meta_json]:
        if not p.exists():
            raise FileNotFoundError(f"Missing expected export output: {p}")

    meta = _load_json(meta_json)
    polar_atoms = list(meta["polar_atoms"])
    n_polar = len(polar_atoms)

    pz = np.load(placements_npz, allow_pickle=True)
    x12 = pz["xform_world_stub_12_f32"]
    irot = pz["irot_id_u16"]
    score = pz["score_f32"]
    sat1 = pz["sat1_i16"]
    sat2 = pz["sat2_i16"]

    if x12.ndim != 2 or x12.shape[1] != 12:
        raise AssertionError(f"Bad xform shape: {x12.shape}")
    n = int(x12.shape[0])
    if not (irot.shape == score.shape == sat1.shape == sat2.shape == (n,)):
        raise AssertionError("Placement arrays must all have shape (N,).")
    if np.any(score > 0.0):
        raise AssertionError("Expected exported scores to be <= 0.0 for rifdock insertion gating.")

    def _sat_ok(v: np.ndarray) -> bool:
        return bool(np.all((v == -1) | ((v >= 0) & (v < n_polar))))

    if not _sat_ok(sat1) or not _sat_ok(sat2):
        raise AssertionError("sat1/sat2 indices out of range for polar_atoms.")

    lz = np.load(lib_npz, allow_pickle=True)
    lib_irot = lz["irot_id_u16"]
    lib_xyz = lz["center_atom_xyz_stub_f32"]
    if lib_irot.ndim != 1:
        raise AssertionError("irot_lib irot_id_u16 must be 1D.")
    if lib_xyz.ndim != 3 or lib_xyz.shape[2] != 3:
        raise AssertionError(f"irot_lib center_atom_xyz_stub_f32 must be (M,n_atoms,3), got {lib_xyz.shape}")

    # All irot ids used in placements must be present in the lib.
    used = set(np.unique(irot).tolist())
    have = set(lib_irot.tolist())
    missing = sorted(used - have)
    if missing:
        raise AssertionError(f"placements reference missing irot ids: {missing[:10]}")


if __name__ == "__main__":
    main()

