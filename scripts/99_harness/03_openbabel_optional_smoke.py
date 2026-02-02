#!/usr/bin/env python
from __future__ import annotations

import argparse
import subprocess
import time
from pathlib import Path


def _write_text(path: Path, s: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(s, encoding="utf-8")


def _run(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True)


def _toy_pdb() -> str:
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
    ap = argparse.ArgumentParser(description="Smoke: scripts must run without importing OpenBabel when not needed.")
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    ap.add_argument("--tag", type=str, default=str(int(time.time())))
    args = ap.parse_args()

    root = args.repo_root.resolve()
    outdir = root / "processed/99_harness" / f"openbabel_optional_{args.tag}"
    outdir.mkdir(parents=True, exist_ok=True)

    ligand = outdir / "LIG.pdb"
    out_sites = outdir / "polar_sites.json"
    _write_text(ligand, _toy_pdb())

    # 1) Polar-site extraction must support backend=rdkit even if OpenBabel is missing/broken.
    _run(
        [
            "uv",
            "run",
            "-p",
            "3.11",
            "python",
            str(root / "scripts/02_polar_sites/01_extract_polar_sites.py"),
            str(ligand),
            "-o",
            str(out_sites),
            "--backend",
            "rdkit",
        ]
    )

    # 2) Candidate generator should at least be invocable for --help without requiring OpenBabel import.
    _run(["uv", "run", "-p", "3.11", "python", str(root / "scripts/04_candidates/01_generate_candidates.py"), "--help"])


if __name__ == "__main__":
    main()
