# Quickstart — vdM × RIFGen deterministic ligand motif

Goal: given a **ligand bound pose** (protein–ligand complex, ligand pose treated as fixed), generate a **deterministic 8–15 residue “motif”** such that **every ligand polar atom is satisfied** (H-bonds / salt bridges) to compensate desolvation loss.

This repo is a lightweight project log + planning space. The actual upstream codebases are vendored under `external/` for code reading:

- `external/Combs2024` @ `c17476366b07e1861a94e133da5bf3a61ba5d6ee`
- `external/rifdock` @ `62c3a1a19438f443a7bc391cce5c0574bea8b02b`

## What you get (deliverables)

1) A **data spec** to convert Combs vdM parquet databases into a compact “vdXform” library that rifgen/rifdock can load efficiently.

2) A **deterministic motif solver spec** (8–15 residues, hard constraint: cover all ligand polar atoms).

3) An **MVP execution path** with minimal rifdock/rifgen changes (Python pre-processing; avoid rifgen randomness).

See `docs/design_plan.md` for the detailed architecture, quantitative constraints, and milestones.
See `docs/vdxform_spec.md` for the current parquet→vdXform conversion spec.

## Required inputs (frozen)

- Protein–ligand complex structure containing the ligand **bound pose**.
- Ligand atom typing (HBD/HBA, cation/anion) + mapping to “chemical-group frames”.

## Local assets you already have

- Ligand chemical-group (CG) detection + atom mapping script:
  - `~/practice_arena/050_test_COMBS/vdm_designer/vdm_core/generate_cg_atommap.py`
  - It detects Combs-style CGs in a ligand (PDB/MOL2) using SMARTS + connectivity ordering and outputs a `vdm_cg_aa_atommap_dict` mapping entries that include:
    - `cg`: CG name (e.g. `coo`, `conh2`, `ccn`, `phenol`, `ph`, `hid/hie/hip`, ...)
    - `lgd_sel`: ordered ligand atom names (must align with `correspond_names`)
    - `correspond_resname` / `correspond_names`: the AA identity + atom order from `cg_dicts.txt`
    - `represent_name`: a representative AA atom (e.g. `OD1`, `NE2`, `NZ`) used as an anchor in some workflows

- Combs2024 vdM database path (on your machine):
  - `~/practice_arena/050_test_COMBS/Combs2024/database/database/`
  - In this dataset: `vdMs/` is ~5.2GB and contains 354 `*.parquet.gzip` files organized by CG.

## Minimal pipeline (MTX example)

### 1) CG mapping (Combs-style)

Example (JSON output used in this repo):

```bash
uv sync -p 3.11 --extra rdkit
uv run -p 3.11 python ~/practice_arena/050_test_COMBS/vdm_designer/vdm_core/generate_cg_atommap.py \
  ligand.pdb \
  --format json \
  --vdm-database ~/practice_arena/050_test_COMBS/Combs2024/database/database/vdMs \
  -o outputs/01_cgmap/ligand_cg_atommap.json
```

Determinism note: this script is deterministic given the same input file and the same bond-perception backend (it may call OpenBabel if available).

### 2) Polar sites + deterministic ligand site frames

For MTX in this repo:

```bash
bash scripts/02_polar_sites/01_run_mtx_polar_sites.sh
```

Outputs:
- `outputs/02_polar_sites/MTX_polar_sites.json` (polar atom list)
- `outputs/02_polar_sites/MTX_site_frames.json` (deterministic iFG frames; includes a `bb_cnh` fallback for uncovered donor atoms)

### 3) Convert Combs parquet → vdXform library (CPU-heavy; run on `dg` if needed)

Debug (small, capped) conversion for the CGs needed by MTX:

```bash
bash scripts/03_vdxform/03_run_mtx_needed_cgs_debug.sh
```

Full conversion (recommended on `dg`):

```bash
UV_WORKERS=48 bash scripts/03_vdxform/02_run_mtx_needed_cgs.sh
```

### 4) Generate residue-placement candidates (MVP)

```bash
bash scripts/04_candidates/01_run_mtx_candidates_debug.sh
```

Outputs:
- `processed/04_candidates/MTX_candidates_debug.npz`
- `processed/04_candidates/MTX_candidates_debug.json`

Non-debug run (expected to be heavier; intended for `dg`):

```bash
TOP_PER_SITE=2000 bash scripts/04_candidates/02_run_mtx_candidates.sh
```

### 5) Solve deterministic 8–15 residue motif (MVP)

```bash
bash scripts/05_solver/01_run_mtx_solver_debug.sh
```

Outputs:
- `outputs/05_solver/MTX_motif_debug.json`
- `outputs/05_solver/MTX_motif_debug.pdb` (full-atom residues where available via `center_atom_xyz_stub_f32`, plus the ligand)

Non-debug run:

```bash
TIME_LIMIT_S=300 bash scripts/05_solver/02_run_mtx_solver.sh
```

### 6) Determinism regression (recommended while iterating)

Runs candidate generation + solver twice and asserts byte-identical outputs:

```bash
bash scripts/99_harness/01_run_mtx_determinism_regression.sh
```

## Environment notes (important)

- This workspace’s system `python3` is 3.6.8, and the system RDKit is broken (missing FreeINCHI symbols).
- Parquet IO for the Combs vdM database requires `pyarrow` or `fastparquet`. In this environment `pyarrow` is not installed, and modern `pyarrow` wheels typically do not support Python 3.6.
- Practical consequence: the **parquet → vdXform** conversion step should run in a dedicated modern Python env installed via **uv** (or pixi).

### uv setup (recommended)

This repo ships a minimal `pyproject.toml` for a modern toolchain (Python ≥ 3.11, `pyarrow`, `ortools`, etc.).

```bash
uv python install 3.11
uv lock -p 3.11
uv sync -p 3.11 --extra rdkit
```

## Core references

- vdM concept (“van der Mer”) — Polizzi & DeGrado, *Science* (2020): “A defined structural unit enables de novo design of small-molecule–binding proteins.” DOI: `10.1126/science.abb8330`.
- RIF / RifGen / RifDock concept (“rotamer interaction field”, hierarchical search) — Cao, Coventry, et al., *Nature* (2022): “Design of protein-binding proteins from the target structure alone.” DOI: `10.1038/s41586-022-04654-9`.
- A RifDock application for small-molecule-enabled protein engineering (example of RIFDock in a small-molecule context) — Langan et al., *Nat. Biotechnol.* (2019): “De novo design of bioactive protein switches.” DOI: `10.1038/s41587-019-0260-5`.

Notes:
- The *Nature* (2022) paper is the most direct “primary literature” description of RifGen/RifDock mechanics (RIF generation as inverse rotamers stored in a 6D hash; branch-and-bound/hierarchical search). It also points to the open-source implementation: `https://github.com/rifdock/rifdock`.
- The *Science* (2020) vdM paper is the primary motivation for using PDB-derived interaction geometries as reusable building blocks.

## How work is tracked

We use **bd (beads)** for multi-session task tracking:

```bash
bd list
bd ready
bd show <id>
```

Key issues created:
- `rifvdm-k8y` (epic)
- `rifvdm-1qi` rifgen scaling/sampling constraints
- `rifvdm-02v` Combs vdM schema + re-encoding decision
- `rifvdm-8np` ligand typing + CG frames
- `rifvdm-50s` vdXform file format + converter
- `rifvdm-kx4` MVP RIF generation path
- `rifvdm-e5b` deterministic motif solver

## Repository layout (current)

- `docs/` project docs (see `docs/changelog.md`, `docs/takeaway.md`)
- `external/` vendored upstream repos for code reading
- `.beads/` bd issue database
