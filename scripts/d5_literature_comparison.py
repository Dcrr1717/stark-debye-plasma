"""Plan D, Phase D5: digit-by-digit comparison against the published
combined-problem literature (the referee's interior-point demand).

Sources (transcribed from the papers, PDFs in references/):
- Yu & Ho, Phys. Plasmas 12, 043302 (2005): Table III/IV (lam_D=10),
  Table VIII (lam_D=5), Table IX (lam_D=inf).  Slater basis, beta=2,
  N=28, lmax=11, CRM.
- Lin & Ho, J. Phys. B 44, 175001 (2011): Table 1 (lam_D=inf, 10,
  ground), Table 2 (lam_D=inf, n=2 doublet at F=0.02, lmax=20).
  Lagrange-Laguerre mesh N=40, complex rotation.

Our map: d3_map.csv (single discretization nb=70, lmax=14, beta=0.35);
n=2: d3_manifold.csv.  lam = 1/lam_D.
"""

import csv
from pathlib import Path

D = Path(__file__).parent

M = {}
for r in csv.DictReader((D / "d3_map.csv").open()):
    M[(round(float(r["F"]), 3), round(float(r["lam"]), 3))] = (
        float(r["Er"]), float(r["Gamma"]) / 2, float(r["qual"]))
C = {}
for r in csv.DictReader((D / "d3_manifold.csv").open()):
    C[(r["comp"], round(float(r["lam"]), 3))] = (
        float(r["Er"]), float(r["Gamma"]) / 2, float(r["qual"]))

# (label, F, lam, Er_lit, halfGamma_lit, source)
GROUND = [
    ("lamD=10 ", 0.02, 0.1, -0.4080089060564, None,      "YuHo T.III"),
    ("lamD=10 ", 0.04, 0.1, -0.4110284061907, 6.376202368342e-6, "YuHo T.III/IV"),
    ("lamD=10 ", 0.04, 0.1, -0.4110284061850, 6.3762009250e-6,   "LinHo T.1"),
    ("lamD=10 ", 0.06, 0.1, -0.4168135614846, 4.324681742989e-4, "YuHo T.III/IV"),
    ("lamD=10 ", 0.08, 0.1, -0.4254044082694, 2.964439958367e-3, "YuHo T.III/IV"),
    ("lamD=10 ", 0.08, 0.1, -0.4254044083057, 2.9644392780e-3,   "LinHo T.1"),
    ("lamD=10 ", 0.10, 0.1, -0.4351823341577, 8.521062774048e-3, "YuHo T.III/IV"),
    ("lamD=5  ", 0.02, 0.2, -0.3278807711360, 3.19780521255e-9,  "YuHo T.VIII"),
    ("lamD=5  ", 0.04, 0.2, -0.3313840162260, 4.65354755588e-5,  "YuHo T.VIII"),
    ("lamD=5  ", 0.06, 0.2, -0.3380688974550, 1.12356514641e-3,  "YuHo T.VIII"),
    ("lamD=5  ", 0.08, 0.2, -0.3470424232710, 5.03126160679e-3,  "YuHo T.VIII"),
    ("lamD=5  ", 0.10, 0.2, -0.3565207954490, 1.19057120686e-2,  "YuHo T.VIII"),
    ("lamD=inf", 0.04, 0.0, -0.5037715910000, 1.94634960e-6,     "YuHo T.IX"),
    ("lamD=inf", 0.08, 0.0, -0.5175606170000, 2.26982913e-3,     "YuHo T.IX"),
    ("lamD=inf", 0.10, 0.0, -0.5274181740000, 7.26905777e-3,     "YuHo T.IX"),
]

print("GROUND STATE vs literature (ours = single-grid map, spread q):")
print(f"{'point':>22} {'dEr':>9} {'dG/2':>9} {'q':>7}  source")
worst_er = worst_g = 0.0
for lab, F, lam, er, hg, src in GROUND:
    ours = M[(F, lam)]
    dEr = abs(ours[0] - er)
    tag = f"{lab} F={F:.2f}"
    dG = abs(ours[1] - hg) if hg is not None else float("nan")
    if hg is not None and hg > 1e-8:
        worst_g = max(worst_g, dG / hg)
    worst_er = max(worst_er, dEr)
    print(f"{tag:>22} {dEr:9.1e} {dG:9.1e} {ours[2]:7.0e}  {src}")
print(f"   worst |dEr| = {worst_er:.1e};  worst rel dG/2 (G/2>1e-8) = {worst_g:.1e}")

# n=2 doublet, lam_D = inf, F=0.02 (Lin-Ho Table 2, lmax=20 values)
print("\nN=2 DOUBLET at F=0.02, lam=0 vs Lin-Ho Table 2 (lmax=20):")
LIT2 = {
    "blue": (-0.0889837856, 0.0154457580),   # upper (k=+1)
    "red": (-0.2066822365, 0.0303928560),    # lower (k=-1)
}
for comp, (er, hg) in LIT2.items():
    ours = C[(comp, 0.0)]
    print(f"   {comp:>4}: ours Er={ours[0]:+.10f} G/2={ours[1]:.10f} (q {ours[2]:.0e})")
    print(f"         lit  Er={er:+.10f} G/2={hg:.10f}"
          f"   dEr={abs(ours[0]-er):.1e}  dG/2={abs(ours[1]-hg):.1e}")

# full-precision our values for the paper table
print("\nOur full-precision values for the comparison table:")
for F, lam in ((0.04, 0.1), (0.08, 0.1), (0.10, 0.1), (0.06, 0.2), (0.10, 0.2)):
    o = M[(F, lam)]
    print(f"   F={F:.2f} lam={lam:.1f}: Er={o[0]:.10f}  G/2={o[1]:.10e}  q={o[2]:.0e}")
for comp in ("red", "blue"):
    o = C[(comp, 0.0)]
    print(f"   n=2 {comp}: Er={o[0]:.10f}  G/2={o[1]:.10e}  q={o[2]:.0e}")
