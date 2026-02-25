#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import statistics
from pathlib import Path
from typing import Any


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _get_param(r: dict[str, Any], keys: list[str], *, cast):
    p = r.get("params", {})
    for k in keys:
        if k in p:
            return cast(p[k])
    raise KeyError(f"Missing parameter in run {r.get('run_id')}: any of {keys}")


def _run_key(r: dict[str, Any]) -> tuple[int, float]:
    min_sc = _get_param(
        r,
        ["min_sidechain_contact_count", "min_pocket_sidechain_contacts"],
        cast=int,
    )
    weight = _get_param(
        r,
        ["sidechain_contact_weight", "pocket_contact_weight"],
        cast=float,
    )
    return (min_sc, weight)


def _build_run_map(summary: dict[str, Any]) -> dict[tuple[int, float], dict[str, Any]]:
    m: dict[tuple[int, float], dict[str, Any]] = {}
    for r in summary.get("runs", []):
        m[_run_key(r)] = r
    return m


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Compare MTX pocket sweep runtimes between two summary.json outputs."
    )
    ap.add_argument("--local-summary", type=Path, required=True)
    ap.add_argument("--dg-summary", type=Path, required=True)
    ap.add_argument("--out-md", type=Path, required=True)
    ap.add_argument("--out-json", type=Path, default=None)
    args = ap.parse_args()

    local = _read_json(args.local_summary)
    dg = _read_json(args.dg_summary)

    local_map = _build_run_map(local)
    dg_map = _build_run_map(dg)
    keys = sorted(set(local_map).intersection(dg_map))
    if not keys:
        raise SystemExit("No overlapping (min_sc_contacts, weight) runs between summaries.")

    rows: list[dict[str, Any]] = []
    local_total_s: list[float] = []
    dg_total_s: list[float] = []

    for key in keys:
        lr = local_map[key]
        dr = dg_map[key]
        lt = float((lr.get("timing_s") or {}).get("total_s", 0.0))
        dt = float((dr.get("timing_s") or {}).get("total_s", 0.0))
        local_total_s.append(lt)
        dg_total_s.append(dt)
        rows.append(
            {
                "min_sc_contacts": int(key[0]),
                "weight": float(key[1]),
                "local_total_s": lt,
                "dg_total_s": dt,
                "speedup_local_over_dg": (lt / dt) if dt > 0 else None,
                "local_candidates": int((lr.get("metrics") or {}).get("n_candidates", 0)),
                "dg_candidates": int((dr.get("metrics") or {}).get("n_candidates", 0)),
            }
        )

    local_mean = statistics.mean(local_total_s)
    dg_mean = statistics.mean(dg_total_s)
    local_median = statistics.median(local_total_s)
    dg_median = statistics.median(dg_total_s)
    mean_speedup = local_mean / dg_mean if dg_mean > 0 else 0.0

    report = {
        "local_summary": str(args.local_summary.resolve()),
        "dg_summary": str(args.dg_summary.resolve()),
        "n_matched_runs": len(rows),
        "rows": rows,
        "aggregate": {
            "local_total_s_mean": local_mean,
            "dg_total_s_mean": dg_mean,
            "local_total_s_median": local_median,
            "dg_total_s_median": dg_median,
            "mean_speedup_local_over_dg": mean_speedup,
        },
    }

    if args.out_json is not None:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    md_lines = [
        "# MTX Pocket Runtime Compare (local vs dg)",
        "",
        f"- local summary: `{Path(report['local_summary'])}`",
        f"- dg summary: `{Path(report['dg_summary'])}`",
        f"- matched runs: {report['n_matched_runs']}",
        "",
        "| min_sc_contacts | weight | local_total_s | dg_total_s | speedup(local/dg) | local_candidates | dg_candidates |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        speed = row["speedup_local_over_dg"]
        speed_s = f"{float(speed):.3f}" if speed is not None else "NA"
        md_lines.append(
            "| "
            + f"{row['min_sc_contacts']} | {row['weight']:.2f} | {row['local_total_s']:.2f} | {row['dg_total_s']:.2f} | "
            + f"{speed_s} | {row['local_candidates']} | {row['dg_candidates']} |"
        )
    md_lines.extend(
        [
            "",
            f"- mean local total_s: {local_mean:.2f}",
            f"- mean dg total_s: {dg_mean:.2f}",
            f"- median local total_s: {local_median:.2f}",
            f"- median dg total_s: {dg_median:.2f}",
            f"- mean speedup(local/dg): {mean_speedup:.3f}x",
        ]
    )

    args.out_md.parent.mkdir(parents=True, exist_ok=True)
    args.out_md.write_text("\n".join(md_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
