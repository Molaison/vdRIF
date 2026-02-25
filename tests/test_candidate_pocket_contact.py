from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest


def _load_candidates_module():
    repo_root = Path(__file__).resolve().parents[1]
    script_path = repo_root / "scripts/04_candidates/01_generate_candidates.py"
    spec = importlib.util.spec_from_file_location("generate_candidates", script_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_sidechain_contact_features_counts_contact_window() -> None:
    module = _load_candidates_module()

    residue_atoms = np.array(
        [
            # candidate 0
            [
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [3.2, 0.0, 0.0],
                [4.6, 0.0, 0.0],
            ],
            # candidate 1
            [
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [5.5, 0.0, 0.0],
                [2.0, 0.0, 0.0],
            ],
        ],
        dtype=np.float64,
    )
    ligand_xyz = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)

    count, min_dist = module._sidechain_contact_features(
        residue_atoms,
        sidechain_idx=[3, 4],
        ligand_xyz=ligand_xyz,
        min_contact_dist=2.8,
        max_contact_dist=4.8,
    )

    assert count.tolist() == [2, 0]
    assert np.allclose(min_dist, np.array([3.2, 2.0], dtype=np.float32))


def test_sidechain_contact_features_handles_empty_sidechain_idx() -> None:
    module = _load_candidates_module()

    residue_atoms = np.zeros((3, 4, 3), dtype=np.float64)
    ligand_xyz = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)

    count, min_dist = module._sidechain_contact_features(
        residue_atoms,
        sidechain_idx=[],
        ligand_xyz=ligand_xyz,
        min_contact_dist=2.8,
        max_contact_dist=4.8,
    )

    assert count.tolist() == [0, 0, 0]
    assert np.isinf(min_dist).all()


def test_sidechain_contact_features_rejects_invalid_distance_window() -> None:
    module = _load_candidates_module()

    residue_atoms = np.zeros((1, 4, 3), dtype=np.float64)
    ligand_xyz = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)

    with pytest.raises(ValueError):
        module._sidechain_contact_features(
            residue_atoms,
            sidechain_idx=[1],
            ligand_xyz=ligand_xyz,
            min_contact_dist=4.0,
            max_contact_dist=4.0,
        )
