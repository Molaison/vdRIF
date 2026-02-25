# openbabel_meta

This directory provides a minimal workspace member for the `openbabel` distribution.

Why it exists:
- The root `pyproject.toml` maps `openbabel` to a uv workspace source.
- `uv.lock` expects an editable package at `vendor/openbabel_meta`.
- Without this directory, `uv sync` fails in clean worktrees before any scripts can run.

This package is metadata-only and depends on:
- `openbabel-wheel==3.1.1.22`

It intentionally does not ship Python modules.
