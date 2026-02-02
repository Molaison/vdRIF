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

## Reality check (2026-01-15)

- The pipeline can generate a **deterministic motif** and validate **ligand–motif** and **motif-internal** vdW overlaps.
- The “simple satisfaction” validator (`scripts/05_solver/03_validate_motif_polar_satisfaction.py`) is **not equivalent to PLIP**, and in this session we repeatedly saw motifs that “passed” simple checks but **failed PLIP**.
- If you care about “perfect interactions”, treat **PLIP** as the ground-truth gate and iterate candidate generation + solver until `scripts/05_solver/06_validate_motif_plip.py` reports `all_satisfied: true`.

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

Note: the default harness uses `processed/03_vdxform/` (debug vdXform). For MTX, that debug library may not contain any viable candidates for the donor atom `N`, causing a hard failure. Use the full libraries instead:

```bash
uv run -p 3.11 python scripts/99_harness/01_mtx_determinism_regression.py \
  --repo-root . \
  --vdxform-dir processed/03_vdxform_full \
  --solver greedy \
  --top-per-site 200 \
  --time-limit-s 30 \
  --tag full_vdx_smoke
```

### 7) Interaction validation (PLIP) — do not skip

Run PLIP on the generated motif PDB and check that every polar atom is contacted by at least one PLIP interaction:

```bash
uv run -p 3.11 python scripts/05_solver/06_validate_motif_plip.py \
  --motif-pdb outputs/05_solver/MTX_motif.pdb \
  --polar-sites outputs/02_polar_sites/MTX_polar_sites.json \
  -o outputs/05_solver/MTX_motif_plip.json \
  --keep-plip-outdir
```

Also always run both clash validators:

```bash
uv run -p 3.11 python scripts/05_solver/04_validate_motif_clashes.py \
  --motif-pdb outputs/05_solver/MTX_motif.pdb -o outputs/05_solver/MTX_motif_ligand_clash.json --ligand-resname MTX

uv run -p 3.11 python scripts/05_solver/05_validate_motif_internal_clashes.py \
  --motif-pdb outputs/05_solver/MTX_motif.pdb -o outputs/05_solver/MTX_motif_internal_clash.json --ligand-resname MTX
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

## Implementation strategy (mental model)

This repo’s MVP is **vdM-driven candidate placement + deterministic set cover**, not a native rifgen run yet:

1) **Ligand typing + CG mapping**: produce a Combs-style CG atom map (`outputs/01_cgmap/*.json`) and a polar-atom list (`outputs/02_polar_sites/*_polar_sites.json`).
2) **Deterministic ligand site frames**: for each CG mapping entry, build an iFG frame on the ligand (`outputs/02_polar_sites/*_site_frames.json`). For symmetric CGs (e.g., `coo` OD1/OD2 and aromatic `ph`/`phenol` CD1/CD2), emit swapped frames to avoid label/handedness artifacts.
3) **vdXform library**: convert Combs parquet vdMs into a compact NPZ containing, per vdM instance, an `X_ifg_to_stub` transform and the **interaction residue** atom coordinates in the rifdock stub frame.
4) **Candidate generation** (`scripts/04_candidates/01_generate_candidates.py`):
   - For each site frame, load the corresponding `vdxform_<cg>.npz`.
   - Compose placement `X_world_stub = X_world_ifg ∘ X_ifg_to_stub`.
   - Fast prefilter by ligand clash with N/CA/C/CB, then compute a **satisfaction bitmask** for the polar atoms assigned to this site.
   - Apply sidechain-facing filtering and a full-atom ligand clash check.
   - Keep top-K per site (deterministic ties by `vdm_id`).
5) **Motif solving** (`scripts/05_solver/01_solve_motif.py`): choose 8–15 candidates that cover all polar bits while avoiding residue–residue clashes; write motif PDB.

## Session takeaways: recurring problems + what we changed

### A) “Looks satisfied” but PLIP says no

Symptom:
- `03_validate_motif_polar_satisfaction.py` reported `all_satisfied: true`, but `06_validate_motif_plip.py` still reported unsatisfied polar atoms.

Root causes we hit:
- **Over-permissive motif atom typing**: we previously treated `MET:SD` and `CYS:SG` as acceptors (and even `GLN:NE2` as acceptor), which can create false-positive “satisfaction” by pure distance.
- **Ligand donor geometry mismatch**: using ligand H coordinates directly from the input PDB can be wrong/unphysical (even if the ligand has explicit Hs), which breaks donor-angle/HA-distance gating and disagrees with PLIP.

What to do:
- Always gate success using `scripts/05_solver/06_validate_motif_plip.py`.
- For PLIP-aligned runs, keep acceptor typing conservative (Met/Cys sulfur are not acceptors in PLIP/OpenBabel): use `--acceptor-model plip` where available.
- For ligand donors, prefer OpenBabel-rebuilt H coordinates to align with PLIP.

### B) Severe “motif-internal” clashes despite passing ligand clash checks

Symptom:
- Motif did not clash with ligand, but residues clashed badly with each other (user-visible “motif not self-compatible”).

Root cause:
- We originally only validated **ligand–motif** clashes. Internal clashes require a separate check.

What to do:
- Always run `scripts/05_solver/05_validate_motif_internal_clashes.py` and fail the run if it reports overlaps.

### C) “Sampling too small” / unexpectedly few candidates

Symptom:
- Only a tiny number of candidates survived, despite vdM libraries being huge.

Root causes we hit:
- Using **debug vdXform** (capped parquet conversion) produces tiny libraries (`processed/03_vdxform/*/vdxform_*.npz` ~1.5MB each). This is only for smoke tests and can make coverage look impossible.
- Aggressive geometric gates + clash gates can wipe out all candidates for some polar atoms (especially ligand donors if H geometry is inconsistent).

What to do:
- Use the full libraries under `processed/03_vdxform_full/` when evaluating feasibility.
- Increase `--top-per-site` and consider `--top-per-site-per-atom` to avoid rare bits getting pruned by a shared heap.

### D) vdM semantics confusion: which residue should become the “motif residue”?

This was a major source of incorrect motifs early on.

Ground truth in Combs2024 parquet (verify in `external/Combs2024/combs2/design/functions.py:get_extended`):
- `chain == 'Y'`: iFG/CG atoms (used for superposition / defines the iFG frame).
- `chain == 'X'`: the **segment in the protein** around the interacting residue (typically resnum 9/10/11).
- The “interaction residue” is **`chain=='X' & resnum==10`** (what we export/place as motif residue in MVP).

If you export/place the iFG residue itself, you can accidentally overlay ligand-matched atoms and generate nonsense.

### E) Scaling mismatch vs rifgen

rifgen assumes a relatively small, fixed rotamer set per amino acid, while vdM libraries can be millions:
- Any approach that “just loops over vdMs” will bottleneck on I/O and on full-atom clash checks.
- The MVP therefore uses (1) per-site top-K truncation, (2) vectorized prefilters, and (3) a deterministic solver.

### F) Beads gotcha: repo mismatch

If `bd ready` warns `DATABASE MISMATCH DETECTED`, run:

```bash
bd migrate --update-repo-id
```

…before syncing, otherwise beads sync can delete/lose issues across clones.
