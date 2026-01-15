#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


def _load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _first_two_bits(mask_u16: int) -> tuple[int, int]:
    bits = []
    for b in range(16):
        if (mask_u16 >> b) & 1:
            bits.append(b)
            if len(bits) == 2:
                break
    if not bits:
        return -1, -1
    if len(bits) == 1:
        return bits[0], -1
    return bits[0], bits[1]


@dataclass(frozen=True)
class IRotKey:
    aa3: str
    cg: str
    cluster_number: int


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Export rifdock-facing inputs from our candidate NPZ: a ligand-specific irot library + placement list.\n"
            "This is an intermediate format so C++/rifdock can build a real .rif using its native XformHash."
        )
    )
    ap.add_argument("--candidates-npz", type=Path, required=True)
    ap.add_argument("--candidates-meta", type=Path, required=True)
    ap.add_argument("--site-frames", type=Path, required=True, help="The ligand site frames JSON used to generate candidates.")
    ap.add_argument("--out-prefix", type=Path, required=True, help="Output prefix (writes *_placements.npz, *_irot_lib.npz, *_meta.json).")
    ap.add_argument(
        "--rotamer-bits",
        type=int,
        default=12,
        help="Max bits for irot ids (default 12 -> 4096). Export fails if unique irots exceed this.",
    )
    ap.add_argument(
        "--score-sign",
        type=str,
        default="negate",
        choices=["negate", "keep"],
        help="rifdock RIF insert ignores score>0; by default we negate our score so exported scores are <=0.",
    )
    args = ap.parse_args()

    meta = _load_json(args.candidates_meta)
    polar_atoms: list[str] = meta["polar_atoms"]

    frames = _load_json(args.site_frames)
    site_cg = [str(s["cg"]) for s in frames["sites"]]

    z = np.load(args.candidates_npz, allow_pickle=True)
    cand_id = z["cand_id_u64"].astype(np.uint64)
    site_index = z["site_index_u16"].astype(np.int64)
    aa3 = z["aa3"].astype("U3")
    score = z["score_f32"].astype(np.float32)
    cover = z["cover_mask_u16"].astype(np.uint16)
    xform12 = z["xform_world_stub_12_f32"].astype(np.float32)
    center_atom_xyz_stub = z["center_atom_xyz_stub_f32"].astype(np.float32)
    cluster_number = z["cluster_number_i32"].astype(np.int32)

    n = int(cand_id.shape[0])
    if n == 0:
        raise ValueError("No candidates provided.")

    cg = np.array([site_cg[int(i)] for i in site_index.tolist()], dtype=object)

    # Build ligand-specific irot ids (deterministic):
    # key = (aa3,cg,cluster_number). This avoids assuming cluster_number is unique across amino acids.
    keys = [IRotKey(str(a), str(c), int(k)) for a, c, k in zip(aa3.tolist(), cg.tolist(), cluster_number.tolist(), strict=True)]
    uniq = sorted({(k.aa3, k.cg, k.cluster_number) for k in keys})

    max_ids = 1 << int(args.rotamer_bits)
    if len(uniq) > max_ids:
        raise RuntimeError(
            f"Too many unique irots for rotamer_bits={args.rotamer_bits}: {len(uniq)} > {max_ids}. "
            "Need stronger compression / different rif_type."
        )

    key_to_irot: dict[tuple[str, str, int], int] = {k: i for i, k in enumerate(uniq)}
    irot_id = np.array([key_to_irot[(k.aa3, k.cg, k.cluster_number)] for k in keys], dtype=np.uint16)

    # Representative coords per irot: choose the best-scoring candidate deterministically.
    best_idx: dict[int, int] = {}
    for i in range(n):
        rid = int(irot_id[i])
        if rid not in best_idx:
            best_idx[rid] = i
            continue
        j = best_idx[rid]
        if float(score[i]) > float(score[j]) or (float(score[i]) == float(score[j]) and int(cand_id[i]) < int(cand_id[j])):
            best_idx[rid] = i

    irot_ids_sorted = np.array(sorted(best_idx.keys()), dtype=np.uint16)
    repr_i = np.array([best_idx[int(r)] for r in irot_ids_sorted], dtype=np.int64)

    # sat1/sat2 from coverage bitmask (best-effort compatibility with RotamerScoreSat).
    sat12 = np.array([_first_two_bits(int(m)) for m in cover.tolist()], dtype=np.int16)
    sat1 = sat12[:, 0]
    sat2 = sat12[:, 1]

    score_out = score.copy()
    if args.score_sign == "negate":
        score_out = -score_out

    out_prefix = args.out_prefix
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    np.savez_compressed(
        out_prefix.with_name(out_prefix.name + "_placements.npz"),
        xform_world_stub_12_f32=xform12,
        irot_id_u16=irot_id,
        score_f32=score_out,
        sat1_i16=sat1.astype(np.int16),
        sat2_i16=sat2.astype(np.int16),
    )

    # irot library: only representatives, but includes the key so C++ can reconstruct per-aa geometry.
    lib_aa3 = np.array([uniq[int(r)][0] for r in irot_ids_sorted.tolist()], dtype="U3")
    lib_cg = np.array([uniq[int(r)][1] for r in irot_ids_sorted.tolist()], dtype=object)
    lib_cluster = np.array([uniq[int(r)][2] for r in irot_ids_sorted.tolist()], dtype=np.int32)
    np.savez_compressed(
        out_prefix.with_name(out_prefix.name + "_irot_lib.npz"),
        irot_id_u16=irot_ids_sorted,
        aa3=lib_aa3,
        cg=lib_cg,
        cluster_number_i32=lib_cluster,
        center_atom_xyz_stub_f32=center_atom_xyz_stub[repr_i],
        repr_cand_id_u64=cand_id[repr_i],
    )

    _write_json(
        out_prefix.with_name(out_prefix.name + "_meta.json"),
        {
            "inputs": {
                "candidates_npz": str(args.candidates_npz),
                "candidates_meta": str(args.candidates_meta),
                "site_frames": str(args.site_frames),
            },
            "params": {
                "rotamer_bits": int(args.rotamer_bits),
                "score_sign": str(args.score_sign),
            },
            "polar_atoms": polar_atoms,
            "sat_index_semantics": {
                "sat1_i16": "index into polar_atoms (or -1)",
                "sat2_i16": "index into polar_atoms (or -1)",
                "note": "Derived from cover_mask_u16 by taking the first two set bits (deterministic).",
            },
            "counts": {
                "n_candidates": int(n),
                "n_irots": int(irot_ids_sorted.shape[0]),
            },
        },
    )


if __name__ == "__main__":
    main()

