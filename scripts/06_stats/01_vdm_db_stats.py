#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import pyarrow.parquet as pq


def _load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


@dataclass
class CgStats:
    cg: str
    n_files: int = 0
    n_row_groups: int = 0
    n_rows: int = 0
    n_instances: int = 0
    n_probe_names: int = 0
    n_clusters: int = 0
    cluster_min: int | None = None
    cluster_max: int | None = None


def _iter_parquet_files(db: Path, cg: str) -> Iterable[Path]:
    d = db / cg
    if not d.exists():
        return []
    return sorted(d.glob("*.parquet.gzip"))


def main() -> None:
    ap = argparse.ArgumentParser(description="Scan Combs2024 vdMs parquet DB and summarize instance/cluster counts.")
    ap.add_argument("--vdm-db", type=Path, required=True, help="Path to Combs2024 database/database/vdMs")
    ap.add_argument("--cgs", type=str, nargs="+", required=True, help="CG directories to scan (e.g., coo ccn conh2 ph)")
    ap.add_argument("-o", "--out", type=Path, required=True)
    ap.add_argument("--max-files-per-cg", type=int, default=0, help="If >0, only scan first N files per CG (debug).")
    ap.add_argument("--batch-rows", type=int, default=200_000, help="Rows per batch for ParquetFile.iter_batches.")
    args = ap.parse_args()

    db = args.vdm_db
    t0 = time.time()

    out: dict[str, Any] = {
        "inputs": {
            "vdm_db": str(db),
            "cgs": list(args.cgs),
            "max_files_per_cg": int(args.max_files_per_cg),
            "batch_rows": int(args.batch_rows),
        },
        "cgs": {},
    }

    for cg in args.cgs:
        st = CgStats(cg=cg)
        files = list(_iter_parquet_files(db, cg))
        if args.max_files_per_cg > 0:
            files = files[: int(args.max_files_per_cg)]

        # These sets are intentionally in-memory: for MTX-relevant CGs we expect O(1e5–1e6) keys.
        instance_keys: set[str] = set()
        probe_names: set[str] = set()
        clusters: set[int] = set()

        for pf in files:
            st.n_files += 1
            pqf = pq.ParquetFile(pf)
            st.n_row_groups += int(pqf.num_row_groups)
            st.n_rows += int(pqf.metadata.num_rows) if pqf.metadata is not None else 0

            want = [c for c in ["rota", "CG", "probe_name", "cluster_number"] if c in pqf.schema.names]
            if "probe_name" not in want or "rota" not in want:
                continue

            for batch in pqf.iter_batches(columns=want, batch_size=int(args.batch_rows)):
                cols = {nm: batch.column(i) for i, nm in enumerate(want)}
                rota = cols["rota"].to_pylist()
                cg_i = cols.get("CG")
                cg_v = cg_i.to_pylist() if cg_i is not None else [None] * len(rota)
                probe = cols["probe_name"].to_pylist()
                cln_i = cols.get("cluster_number")
                cln_v = cln_i.to_pylist() if cln_i is not None else [None] * len(rota)

                # Group key in Combs is (rota, CG, probe_name). Using a compact string key is fast enough.
                for r, cgi, pn, cln in zip(rota, cg_v, probe, cln_v):
                    if pn is not None:
                        probe_names.add(str(pn))
                    instance_keys.add(f"{int(r)}|{int(cgi) if cgi is not None else -1}|{pn}")
                    if cln is None:
                        continue
                    try:
                        cln_int = int(cln)
                    except Exception:
                        continue
                    clusters.add(cln_int)

        st.n_instances = len(instance_keys)
        st.n_probe_names = len(probe_names)
        st.n_clusters = len(clusters)
        if clusters:
            st.cluster_min = min(clusters)
            st.cluster_max = max(clusters)

        out["cgs"][cg] = {
            "cg": st.cg,
            "n_files": st.n_files,
            "n_row_groups": st.n_row_groups,
            "n_rows": st.n_rows,
            "n_instances": st.n_instances,
            "n_probe_names": st.n_probe_names,
            "n_clusters": st.n_clusters,
            "cluster_min": st.cluster_min,
            "cluster_max": st.cluster_max,
        }

    out["runtime_s"] = float(time.time() - t0)
    _write_json(args.out, out)


if __name__ == "__main__":
    main()

