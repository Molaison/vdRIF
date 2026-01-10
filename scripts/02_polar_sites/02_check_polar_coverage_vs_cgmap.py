#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def _load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Check whether cg_atommap covers all RDKit-identified ligand polar atoms."
    )
    ap.add_argument("--polar-sites", type=Path, required=True)
    ap.add_argument("--cg-atommap", type=Path, required=True)
    ap.add_argument("-o", "--out", type=Path, required=True)
    args = ap.parse_args()

    polar = _load_json(args.polar_sites)
    cgmap = _load_json(args.cg_atommap)

    polar_sites = polar["sites"]
    polar_names = {s["atom_name"] for s in polar_sites}

    mapped_names: set[str] = set()
    # cgmap is typically a dict: { "i_j": { "cg": ..., "lgd_sel": [...], ... }, ... }
    if isinstance(cgmap, dict):
        entries = cgmap.values()
    elif isinstance(cgmap, list):
        entries = cgmap
    else:
        raise TypeError(f"Unexpected cg_atommap JSON type: {type(cgmap)}")

    for entry in entries:
        if not isinstance(entry, dict):
            continue
        for name in entry.get("lgd_sel", []) or []:
            mapped_names.add(name)

    covered = sorted(polar_names & mapped_names)
    uncovered = sorted(polar_names - mapped_names)

    payload = {
        "inputs": {"polar_sites": str(args.polar_sites), "cg_atommap": str(args.cg_atommap)},
        "n_polar_atoms": len(polar_names),
        "n_mapped_atom_names": len(mapped_names),
        "n_covered": len(covered),
        "n_uncovered": len(uncovered),
        "covered_atom_names": covered,
        "uncovered_atom_names": uncovered,
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
