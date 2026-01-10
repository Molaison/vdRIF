#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import os
import tempfile
from pathlib import Path
from typing import Any

import numpy as np


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser(description="Pack vdXform parts (*.npz) into a single concatenated npz.")
    ap.add_argument("--parts-dir", type=Path, required=True)
    ap.add_argument("-o", "--out", type=Path, required=True)
    ap.add_argument("--resume", action="store_true")
    args = ap.parse_args()

    if args.resume and args.out.exists():
        return

    part_files = sorted(args.parts_dir.glob("*.npz"))
    if not part_files:
        raise FileNotFoundError(f"No parts found in {args.parts_dir}")

    arrays: dict[str, list[np.ndarray]] = {}
    n_total = 0
    for pf in part_files:
        data = np.load(pf, allow_pickle=True)
        n = int(data["vdm_id_u64"].shape[0])
        n_total += n
        for k in data.files:
            arrays.setdefault(k, []).append(data[k])

    packed: dict[str, np.ndarray] = {}
    for k, parts in arrays.items():
        if k == "aa":
            packed[k] = np.concatenate(parts, axis=0).astype(object)
        else:
            packed[k] = np.concatenate(parts, axis=0)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        dir=str(args.out.parent),
        prefix=args.out.name + ".",
        suffix=".tmp.npz",
        delete=False,
    ) as tf:
        tmp_path = Path(tf.name)
    try:
        np.savez_compressed(tmp_path, **packed)
        os.replace(tmp_path, args.out)
    finally:
        if tmp_path.exists():
            tmp_path.unlink(missing_ok=True)

    _write_json(
        args.out.with_suffix(".json"),
        {"parts_dir": str(args.parts_dir), "n_parts": len(part_files), "n_entries": n_total, "out": str(args.out)},
    )


if __name__ == "__main__":
    main()
