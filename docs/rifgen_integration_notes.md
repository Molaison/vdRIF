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
- `Rot12ScoreSat96`: `RotamerScores<14, RotamerScoreSat<uint32_t,12,...>>` (4096 rotamer ids; larger bins; requires patch `external/patches/rifdock/0001-add-Rot12ScoreSat96.patch`)
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

## Concrete interface proposal (where code changes land)

This repo’s current MVP keeps Python as the correctness oracle and treats rifdock as a **fast, native 6D hash + top‑K container**. That suggests two integration tiers:

### Tier 0 (current MVP, no rifdock C++ changes)

- Python exports an intermediate “RIF input” bundle:
  - placement list: `(X_world_stub, irot_id, score, sat1/sat2)`
  - irot library sidecar: `irot_id -> (aa3,cg,cluster_number, stub-local atoms)`
- See `docs/vdxform_to_rif_spec.md` and `scripts/06_rif_export/01_export_rif_inputs.py`.

This tier avoids reimplementing rifdock’s hasher and packed storage in Python, but it is **not** a `.rif.gz` yet.

### Tier 1 (recommended next step): add a “placement importer” that writes a real `.rif.gz`

Goal: convert the intermediate placement bundle into a rifdock `RifBase` by using rifdock’s own `XformHash` and per-voxel `RotamerScores` insertion.

Where:
- App entry: `external/rifdock/apps/rosetta/rifgen.cc`
  - This already wires `RifFactory` + `RifAccumulator` and writes `.rif.gz`.
  - It currently generates placements via built-in generators (hbonds/apohotspots).
- Insert path: `external/rifdock/apps/rosetta/riflib/rif/RifAccumulators.hh`
  - `RIFAccumulatorMapThreaded::insert(xform, score, rot, sat1, sat2, ...)`
  - Uses `xmap_ptr_->hasher_.get_key(x)` and `value.add_rotamer(rot, score, sat1, sat2, ...)`.

Minimal patch shape:
- Add a new mode to `rifgen.cc` (or a new small app beside it) that:
  1) loads `*_placements.npz` (xform12 + irot_id + score + sat1/sat2)
  2) constructs `rif_accum = rif_factory->create_rif_accumulator(hash_cart_resl, hash_angle_resl, cart_bound, ...)`
  3) loops placements and calls `rif_accum->insert(x, score, irot_id, sat1, sat2, force=false, single_thread=...)`
  4) `rif_accum->condense()` and write `.rif.gz`

Notes:
- Insertion is **hard-gated** by `score <= 0` (see `RIFAccumulatorMapThreaded::insert`: `if( score > 0.0 ) return;`), so exported scores must be non-positive (Python export defaults to negating).
- `irot_id` is stored as an integer field; rifdock does not interpret it during writing. Interpretation happens later when *using* the RIF (rotamer library / atom coords). For vdM motifs we plan to keep `irot_lib` as the authoritative sidecar until/unless we build a custom runtime that can realize these “vdM rotamers”.

### Tier 2 (future): native `RifGeneratorVdM`

If we want to skip Python enumeration and generate directly inside rifdock, a new generator would implement:
- `external/rifdock/apps/rosetta/riflib/rif/RifGenerator.hh` (`RifGenerator::generate_rif(...)`)

It would need:
- a vdM/irot library reader (not parquet; precompiled binary/NPZ->custom C++ format)
- a satisfaction model consistent with current Python oracle
- deterministic chunking/resume and stable tie-breaks

## Where discretization and rotamer-library selection live

### Discretization knobs

- Options and grid ladder live in: `external/rifdock/apps/rosetta/rifgen.cc`
  - `-rifgen:hash_cart_resls`, `-rifgen:hash_ang_resls`, `-rifgen:hash_cart_bounds`
  - `rif_factory->create_rif_accumulator(hash_cart_resl, hash_angle_resl, cart_bound, ...)`

### Rotamer library (`RotamerIndex`) construction

For built-in generators, rifgen builds a `RotamerIndex` from a `RotamerIndexSpec`:
- Spec type: `external/rifdock/schemelib/scheme/chemical/RotamerIndex.hh` (`RotamerIndexSpec`)
- Defaults and cache loading: `external/rifdock/apps/rosetta/riflib/RotamerGenerator.hh` / `.cc`
  - `get_rotamer_spec_default(rot_index_spec, extra_rotamers, ...)`
  - `get_rotamer_index(rot_index_spec, build_per_thread_rotamers)`
  - or `get_rotamer_index(cachefile, build_per_thread_rotamers, rot_index_spec)` to load a cached spec

Some generators can further tweak the spec:
- Example: `external/rifdock/apps/rosetta/riflib/rif/RifGeneratorUserHotspots.cc` (`modify_rotamer_spec(...)`)

For Tier 1 (placement importer), we can bypass `RotamerIndex` entirely because we are not *scoring* placements, only storing them into the RIF container.
