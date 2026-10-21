# Paper D — The Stark problem in a Debye plasma: a coupled-channel treatment

Reproducibility repository for the Paper D manuscript (C. D. Reinoso, 2026),
prepared for submission (Stark resonances under Debye screening).

## Contents
- `main.tex`, `refs.bib` — manuscript source (LaTeX, `iopjournal` class; the
  required `iopjournal.cls` is included).
- `figD_map.pdf`, `figD_manifold.pdf` — figures.
- `data/` — production grid, n=2 manifold, θ trajectories, field-halving check,
  anchors and literature-comparison tables, the m=1 probe data, and
  `settings_production.json`.
- `reproducibility/` — reproducibility notes (`README.md`).
- `scripts/` — the production scripts: solver architecture
  (`d0_architecture.py`), production map (`d1_d3_production.py`,
  `d3_map_v2.py`), referee supplements (`d4_referee_supplements.py`) and
  literature comparison (`d5_literature_comparison.py`).

## How to reproduce
1. Create a Python environment (NumPy/SciPy/matplotlib).
2. Install the spectral engine (see below) so that `scripts/` can import it.
3. Run the scripts in `scripts/` to regenerate the `data/` tables and figures;
   rebuild the PDF with `latexmk -pdf main.tex`.

## Spectral engine
The scripts in `scripts/` depend on the **spectral engine**, which is
maintained as a **separate repository** and must be **cited by its archived
DOI**. Install it into your environment (or place it on `PYTHONPATH`) before
running the scripts; it is deliberately **not** vendored into this repository.
