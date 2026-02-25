#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _run(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True)


def _pack_identity_xforms(translations: np.ndarray) -> np.ndarray:
    n = int(translations.shape[0])
    x12 = np.zeros((n, 12), dtype=np.float32)
    x12[:, 0] = 1.0
    x12[:, 5] = 1.0
    x12[:, 10] = 1.0
    x12[:, 3] = translations[:, 0]
    x12[:, 7] = translations[:, 1]
    x12[:, 11] = translations[:, 2]
    return x12


def _toy_ligand_pdb() -> str:
    return "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n"


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Smoke: balanced objective should produce fuller/more diverse pockets than compact objective."
    )
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    ap.add_argument("--tag", type=str, default=str(int(time.time())))
    args = ap.parse_args()

    root = args.repo_root.resolve()
    outdir = root / "processed/99_harness" / f"solver_balanced_objective_{args.tag}"
    outdir.mkdir(parents=True, exist_ok=True)

    ligand_pdb = outdir / "LIG.pdb"
    candidates_npz = outdir / "toy_candidates.npz"
    candidates_meta = outdir / "toy_candidates.json"
    compact_json = outdir / "compact_report.json"
    compact_pdb = outdir / "compact_motif.pdb"
    balanced_json = outdir / "balanced_report.json"
    balanced_pdb = outdir / "balanced_motif.pdb"
    summary_json = outdir / "summary.json"

    ligand_pdb.write_text(_toy_ligand_pdb(), encoding="utf-8")

    # Synthetic candidate set:
    # - site 0 has high-score broad-coverage options
    # - sites 1..4 are lower-score single-bit options
    # Compact objective will pick min residues, typically concentrated at site 0.
    # Balanced objective should choose more residues and more unique sites.
    polar_atoms = ["A1", "A2", "A3", "A4"]
    atom_order = ["N", "CA", "C", "CB"]
    _write_json(candidates_meta, {"polar_atoms": polar_atoms, "atom_order": atom_order})

    n = 6
    cand_id = np.arange(1, n + 1, dtype=np.uint64)
    site_index = np.array([0, 0, 1, 2, 3, 4], dtype=np.uint16)
    cover_mask = np.array([0b1111, 0b1111, 0b0001, 0b0010, 0b0100, 0b1000], dtype=np.uint16)
    score_f32 = np.array([12.0, 11.0, 8.0, 8.0, 8.0, 8.0], dtype=np.float32)
    aa3 = np.array(["LYS"] * n, dtype="U3")

    # Place candidates far apart to avoid residue-residue clashes in this toy setup.
    translations = np.array(
        [
            [0.0, 0.0, 0.0],
            [20.0, 0.0, 0.0],
            [0.0, 20.0, 0.0],
            [0.0, 0.0, 20.0],
            [20.0, 20.0, 0.0],
            [0.0, 20.0, 20.0],
        ],
        dtype=np.float32,
    )
    xform12 = _pack_identity_xforms(translations)

    center_atom_local = np.array(
        [
            [2.80144, -0.992889, -1.52486],  # N
            [1.95280, 0.220007, -1.52486],  # CA
            [2.87767, 1.43290, -1.52486],  # C
            [1.0264273, 0.25245885, -0.308907],  # CB
        ],
        dtype=np.float32,
    )
    center_atom_xyz_stub = np.broadcast_to(center_atom_local[None, :, :], (n, 4, 3)).copy()

    np.savez_compressed(
        candidates_npz,
        cand_id_u64=cand_id,
        site_index_u16=site_index,
        vdm_id_u64=np.arange(100, 100 + n, dtype=np.uint64),
        cluster_number_i32=np.arange(n, dtype=np.int32),
        aa3=aa3,
        score_f32=score_f32,
        cover_mask_u16=cover_mask,
        xform_world_stub_12_f32=xform12,
        center_atom_xyz_stub_f32=center_atom_xyz_stub,
    )

    solver_py = root / "scripts/05_solver/01_solve_motif.py"
    common = [
        str(solver_py),
        "--candidates-npz",
        str(candidates_npz),
        "--candidates-meta",
        str(candidates_meta),
        "--ligand-pdb",
        str(ligand_pdb),
        "--min-res",
        "2",
        "--max-res",
        "6",
        "--time-limit-s",
        "10",
        "--num-workers",
        "1",
        "--grid-size",
        "4.0",
        "--ca-prefilter",
        "6.0",
        "--clash-tol",
        "0.5",
    ]

    _run(
        [
            sys.executable,
            *common,
            "--out-json",
            str(compact_json),
            "--out-pdb",
            str(compact_pdb),
            "--objective-mode",
            "compact",
        ]
    )

    _run(
        [
            sys.executable,
            *common,
            "--out-json",
            str(balanced_json),
            "--out-pdb",
            str(balanced_pdb),
            "--objective-mode",
            "balanced",
            "--target-res",
            "4",
            "--target-res-penalty",
            "2000",
            "--site-diversity-reward",
            "4000",
        ]
    )

    compact = _load_json(compact_json)
    balanced = _load_json(balanced_json)

    if not bool(compact.get("coverage_complete", False)):
        raise AssertionError("Compact objective failed coverage in toy smoke.")
    if not bool(balanced.get("coverage_complete", False)):
        raise AssertionError("Balanced objective failed coverage in toy smoke.")
    if int(compact["n_selected"]) != 2:
        raise AssertionError(f"Expected compact objective to select 2 residues; got {compact['n_selected']}.")
    if int(compact["n_unique_sites_selected"]) != 1:
        raise AssertionError(
            f"Expected compact objective to use one site in toy setup; got {compact['n_unique_sites_selected']}."
        )
    if int(balanced["n_selected"]) < 4:
        raise AssertionError(f"Expected balanced objective to select at least 4 residues; got {balanced['n_selected']}.")
    if int(balanced["n_unique_sites_selected"]) <= int(compact["n_unique_sites_selected"]):
        raise AssertionError(
            "Balanced objective did not improve unique site usage "
            f"({balanced['n_unique_sites_selected']} <= {compact['n_unique_sites_selected']})."
        )

    _write_json(
        summary_json,
        {
            "compact": {
                "n_selected": int(compact["n_selected"]),
                "n_unique_sites_selected": int(compact["n_unique_sites_selected"]),
                "selected_site_indices": compact.get("selected_site_indices"),
            },
            "balanced": {
                "n_selected": int(balanced["n_selected"]),
                "n_unique_sites_selected": int(balanced["n_unique_sites_selected"]),
                "selected_site_indices": balanced.get("selected_site_indices"),
            },
        },
    )


if __name__ == "__main__":
    main()
