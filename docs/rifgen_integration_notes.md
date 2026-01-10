# rifgen integration notes (vdM-driven ligand motifs)

This note is a pointer map for how rifdock/rifgen works internally, and how we can reuse the same *core idea* (6D discretization + per-bin top‑K) for **vdM‑driven** ligand motif generation.

## rifgen mental model (what it computes)

rifgen builds a **RIF (Rotamer Interaction Field)**: a 6D hash map from a **residue backbone “stub” transform** (translation + orientation) to a small fixed list of **top-scoring rotamers** at that pose.

Key knobs:
- `hash_cart_resl`: translation bin size (Å)
- `hash_angle_resl`: orientation bin size (deg)
- `cart_bound`: spatial bound around the target

Stored per hash bin:
- `K` best entries (`RotamerScores<K,...>`) with (rotamer_id, score, and optional satisfaction bookkeeping).

## Where `RotamerBits` and top‑K come from

rifdock’s default RIF types are instantiated in:
- `external/rifdock/apps/rosetta/riflib/RifFactory.cc`

The packed entry type is `RotamerScore` / `RotamerScoreSat` in:
- `external/rifdock/schemelib/scheme/objective/storage/RotamerScores.hh`

Important: `RotamerBits` is a **compile-time cap** on the maximum number of distinct `rotamer_id` values:
- 9 bits ⇒ 512 max
- 10 bits ⇒ 1024 max

Selected `rif_type` variants (from `RifFactory.cc`):
- `RotScore64`: `RotamerScores<28, RotamerScore<uint16_t,9,...>>`
- `RotScoreSat`: `RotamerScores<14, RotamerScoreSat<uint16_t,9,...>>`
- `Rot10Score6Sat16`: `RotamerScores<14, RotamerScoreSat<uint16_t,10,...>>`
- `Rot12ScoreSat96`: `RotamerScores<14, RotamerScoreSat<uint32_t,12,...>>` (4096 rotamer ids; larger bins)
- `RotScoreSat_1x16`: `RotamerScores<14, RotamerScoreSat<uint16_t,9,..., SatDatum=uint16_t, NSat=1>>`
- `RotScoreSat_2x16`: `RotamerScores<19, RotamerScoreSat<uint16_t,9,..., SatDatum=uint16_t, NSat=2>>`

## What “satisfaction” means inside rifdock

`RotamerScoreSat` stores **up to NSat satisfied target group indices** (`sat1`, `sat2`, ...), not a bitmask.
This is how rifdock enforces `require_satisfaction` constraints for hbond-like interactions.

The scoring path that computes sats lives in:
- `external/rifdock/apps/rosetta/riflib/ScoreRotamerVsTarget.hh` (`score_rotamer_v_target_sat(...)`)

## Non-determinism note (UserHotspots)

`RifGeneratorUserHotspots` samples random perturbations by default (seeded by time), so it is non-deterministic unless you clamp the DOFs:
- `hotspot_sample_cart_bound=0`
- `hotspot_sample_angle_bound=0`
- `hotspot_nsamples=1`

## Mapping to vdM motifs (proposed)

We want the same outer structure:
1) discretize 6D stub poses around the ligand
2) store top‑K **vdM-derived** “rotamers” per bin

But we must replace/extend these assumptions:
- rifgen’s “rotamers” are a small fixed library; vdMs are **10^5–10^6+** instances.
- rifdock packs `rotamer_id` into 9–10 bits; Combs `cluster_number` is **1k–4k** for MTX-relevant CGs.

Therefore, the integration hinges on a bounded **irot library**:
- irot id = a compact identifier for a vdM cluster/representative (≤512/1024)
- store per irot: residue atoms in stub-local coords + vdM prior score(s) + “which ligand polar atom(s) can be satisfied” as sat indices

See `docs/irot_compression_plan.md` and `docs/vdxform_to_rif_spec.md`.
