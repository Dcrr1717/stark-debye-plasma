# Paper D — reproducibility notes

The solver and all generators live in the separate repository
`https://github.com/Dcrr1717/feast-spectral`. The manuscript folder only
carries the exported data (`../data/`).

## Scripts (repo paths, all under `exploration/stark2d/`)

| Script | Role |
|---|---|
| `d0_architecture.py` | `StarkChannels` class: block-tridiagonal coupled-channel assembly on the exact orthonormal Laguerre basis (`src/basis.py`), analytic complex scaling per block, overlap-identity certification. |
| `d1_d3_production.py` | Pure-Stark anchors (two discretizations, two angles) and screening-limit anchor. |
| `d3_map_v2.py` | The 85-point `(F, λ)` ground-state map (→ `data/grid_production.csv`) and the n=2 doublet tracking (→ `data/manifold_n2_F0_02.csv`). Continuation predictor = exact F=0 Yukawa level + previous complex field shift; acceptance = two-angle θc stability (0.18/0.28, retry 0.30/0.42). |
| `d4_referee_supplements.py` | (S1) full θc-trajectories over [0.15, 0.45] for the two least favorable points (→ `data/theta_trajectories.csv`); (S2) quadrature stress test at λ=0.8, θc=0.28 (3000 vs 6000 nodes); (S3) field-shift scaling at λ=0.8, including the F=0.01 halving point (→ `data/field_halving_F001.csv`); (S4) sub-barrier tunneling-exponent check. |
| `d5_literature_comparison.py` | Yu–Ho / Lin–Ho comparison table (→ `data/table_lit_comparison.csv`). 15 ground-state comparisons at λ_D = 10, 5, ∞ giving \|ΔE_r\| ≤ 1.9e-9 and \|ΔΓ/2\| ≤ 3e-9, plus the n=2 doublet against Lin–Ho at matched ℓmax = 14 (all 10 printed digits). |
| `d6_m1_probe.py` | m = 1 probe (→ `data/probe_m1.csv`), with the R1/R2/C1/C2 guards described below. Adds an **optional** `m` argument to `StarkChannels`; m = 0 remains the default, and R2 enforces that the published m = 0 numbers are unchanged bit-for-bit. |

Certification protocol: every constructed basis is checked against the
overlap identity (worst case 1.3e-12); the hydrogen limit of the assembled
system is checked against −1/(2n²) (deviation < 1e-10) before every
production run, and a failed check aborts the run.

## θc-trajectory supplement (editorial task D-S4) — EXPORTED

The trajectories for the two least favorable points (broad ground state at
F=0.02, λ=0.80; positive-energy n=2 component at λ=0.30) are computed by
`d4_referee_supplements.py` (section S1) and are now **stored as a tabular
file**, `../data/theta_trajectories.csv` (14 rows: 2 cases × 7 angles,
θc = 0.15 … 0.45 in steps of 0.05; columns `case, theta_c, Er, Gamma`).

The summary numbers quoted in the manuscript are reproducible from that
file directly:

| Manuscript statement | Recomputed from the CSV |
|---|---|
| ground state F=0.02, λ=0.80 "stationary to 2e-9" | max − min of `Er` over θc ∈ [0.15, 0.45] = **2.21e-9** |
| n=2 uphill λ=0.30 "residual variation 4e-5" | max − min of `Er` over the full range = **3.76e-5** |
| n=2 uphill λ=0.30 "true plateau for θc ≥ 0.25" | max − min of `Er` restricted to θc ≥ 0.25 = **1.15e-8**, a collapse of more than three orders relative to the full range |

No curves were digitized and no values were fabricated: the CSV is the
solver's own S1 output.

## Field-halving check (F = 0.01) — EXPORTED

Section S3 of the same script computes the single off-grid point that backs
the manuscript's non-perturbative-susceptibility argument at λ = 0.8. It is
exported as `../data/field_halving_F001.csv`. The manuscript's factor is
the ratio of the two `shift` columns, 0.0061934258 / 0.0027639468 =
**2.2408** (quoted as 2.24; quadratic response would give 4). The F = 0.02
row of that file reproduces the corresponding `grid_production.csv` row to
all digits.

## m = 1 probe (editorial task D-N1) — NEW PHYSICS, AUTHOR SIGN-OFF REQUIRED

`d6_m1_probe.py` generalizes `StarkChannels` with an **optional** `m`
(default 0, so the published m = 0 path is untouched) and runs a
seven-point probe of the 2p, |m| = 1 resonance at F = 0.02. Output:
`../data/probe_m1.csv`. The script runs four regression/certification
guards *before* the probe and aborts if the two hard ones fail. Verdicts
from the run of record (re-run end to end; the CSV is bit-for-bit
reproducible):

| Guard | What it checks | Verdict |
|---|---|---|
| **R1** | the generalised Gaunt factor `c_l(m) = sqrt(((l+1)²−m²)/((2l+1)(2l+3)))` reduces to the published `(l+1)/sqrt((2l+1)(2l+3))` at m = 0 | **PASS** — the factored form used in the code is **20/20 bit-identical**, max abs diff 0.0e+00. (Recorded for honesty: the *unfactored* expression is only 14/20 bit-identical, max diff 1.1e-16 = 1 ulp. It is **not** the form used.) |
| **R2** | a published production point of `d3_map.csv` (F = 0.02, λ = 0.80, θc = 0.18/0.28) is reproduced **bit-for-bit** after the patch; the script raises `SystemExit` if not | **PASS — BIT-IDENTICAL.** dEr = +0.000e+00, dGamma = +0.000e+00, dspread = +0.000e+00 (Er = −0.050897730206523634, Gamma = 0.013445358359365252, spread = 1.5285630517543646e-10, all three matching the published values exactly). This is the guard that matters: the m-generalisation did not perturb any published number. |
| **C1** | the m = 1 sector reproduces −1/(2n²) with the correct \|m\| = 1 degeneracies (n − 1 states per n) | **PASS** — lowest 6 levels −0.125000000, −0.055555556 (×2), −0.031250000 (×3); max deviation **1.6e-14** |
| **C2** | the m = 1 linear Stark splitting of the n = 3 manifold reproduces the exact 9F — the check that actually exercises `c_l(m=1)` | **PASS** — splitting/F = 9.001051 at F = 1e-4 (rel dev 1.2e-4) and 9.000010 at F = 1e-5 (rel dev 1.2e-6), converging on the exact 9 as F → 0 |

Basis for the probe: ℓmax = 14 → channels ℓ = 1…14, nb = 70, DOF = 980,
overlap certification 7e-13. All seven points were accepted at the
production angle pair 0.18/0.28 (no retry at 0.30/0.42 was needed), with
two-angle spreads from 5e-9 to 2e-6.

**Limits of the probe, stated plainly.** One field row (F = 0.02), one
discretization, seven λ points. The F = 0 Yukawa 2p anchor unbinds near
λ ≈ 0.22, so the last two points (λ = 0.25, 0.30) are placed by linear
extrapolation of the previous two rather than by the Yukawa anchor, and are
correspondingly the least anchored. Most importantly, **there is no
independent literature value for the m = 1 Stark resonance under Debye
screening anywhere in this paper's bibliography**, so C1 and C2 (exact
hydrogenic limits) are the only external anchors that exist for it. The
probe certifies that the m = 1 sector is *correctly assembled*; it does not
claim a converged m = 1 map.

`NEEDS_HUMAN`: the m = 1 probe is new physics added under written author
authorisation. It requires the author's sign-off before submission, and the
decision of whether to keep it in the manuscript, move it to supplementary
material, or drop it, is the author's.

`NEEDS_HUMAN: DOI` — archive scripts + data with a persistent identifier
before submission.
