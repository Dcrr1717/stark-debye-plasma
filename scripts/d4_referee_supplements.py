"""Plan D, Phase D4: referee-demanded supplements.

S1  theta-trajectory at the two least-trustworthy points (broad ground
    at F=0.02, lam=0.8; above-threshold blue n=2 at lam=0.30): show the
    plateau explicitly instead of relying on two angles.
S2  Quadrature stress test for the oscillatory screened blocks
    e^{-lam r e^{i theta}} at the largest lam*theta used: production
    grid (nq=3000) vs doubled grid (nq=6000), max element difference of
    the l=0 screened-potential block and the dipole block.
S3  Quadratic-scaling test of the field shift at lam=0.8: ground at
    F=0.01 vs F=0.02 (ratio 4 = quadratic regime -> polarizability
    language defensible at weak field).
S4  Tunneling-exponent comparison restricted to its sub-barrier domain
    (F < Er^2/4): at F=0.02 that is lam <= 0.30.
"""

import csv
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "exploration" / "stark2d"))

from d0_architecture import StarkChannels  # noqa: E402

t0 = time.perf_counter()
S = StarkChannels(70, 14, 0.35)
print(f"[production basis, cert={S.cert:.0e}]")


def nearest(F, lam, theta, target):
    ev = np.linalg.eigvals(S.build(1.0, lam, F, theta))
    return ev[np.argmin(np.abs(ev - target))]


# ---------------- S1: theta trajectories --------------------------------
print("\nS1 theta-trajectories (production basis):")
cases = [
    ("ground F=0.02 lam=0.80", 0.02, 0.80, -0.0508977 - 0.0067223j),
    ("blue n=2 F=0.02 lam=0.30", 0.02, 0.30, +0.0836898 - 0.0349789j),
]
S1_ROWS = []                     # (D-S4) tabular export of the S1 block
for name, F, lam, seed in cases:
    traj = []
    e = seed
    for th in (0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45):
        e = nearest(F, lam, th, e)
        traj.append((th, e))
        S1_ROWS.append(dict(case=name, theta_c=th,
                            Er=e.real, Gamma=-2 * e.imag))
        print(f"   {name}  th={th:.2f}: {e.real:+.8f}{e.imag:+.8f}j")
    es = [e for _, e in traj]
    # plateau = the 0.18-0.30 production window plus margins
    core = [e for th, e in traj if 0.15 <= th <= 0.35]
    drift = max(abs(a - b) for a in core for b in core)
    print(f"   -> max pairwise drift over th in [0.15,0.35]: {drift:.1e}\n")

# (D-S4) export the S1 trajectories; Gamma = -2 Im E (full width), the
# same convention as d3_map.csv.  Purely additive: no value above is
# touched, only written out at full round-trip precision.
with (Path(__file__).with_name("d4_theta_trajectories.csv")).open(
        "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["case", "theta_c", "Er", "Gamma"])
    w.writeheader()
    for r in S1_ROWS:
        w.writerow(r)

# ---------------- S2: quadrature stress test ----------------------------
print("S2 screened-block quadrature: nq=3000 vs nq=6000 (lam=0.8, th=0.28):")
S2a = StarkChannels(70, 1, 0.35)               # nq = 3000 (default)
S2b = StarkChannels(70, 1, 0.35, nquad=6000)
lam, th = 0.8, 0.28
b = 0.35
for tag, SS in (("nq=3000", S2a), ("nq=6000", S2b)):
    pass
blocks = {}
for tag, SS in (("a", S2a), ("b", S2b)):
    r = SS.x / (2 * b)
    f = np.exp(-lam * r * np.exp(1j * th))
    qw = SS.q[0] * SS.w
    V = -2 * b * ((qw * f / SS.x) @ SS.q[0].T)
    R = ((SS.q[0] * SS.w * r) @ SS.q[1].T)
    blocks[tag] = (V, R)
dV = np.max(np.abs(blocks["a"][0] - blocks["b"][0]))
dR = np.max(np.abs(blocks["a"][1] - blocks["b"][1]))
print(f"   max |dV| screened l=0 block: {dV:.1e}")
print(f"   max |dR| dipole block:       {dR:.1e}")

# ---------------- S3: quadratic scaling of the shift at lam=0.8 ---------
print("\nS3 field-shift scaling at lam=0.8:")
eyuk = np.linalg.eigvalsh(0.5 * (S.build(1.0, 0.8, 0.0, 0.0).real
                                 + S.build(1.0, 0.8, 0.0, 0.0).real.T)).min()
sh = {}
S3_ROWS = []                     # (D-S4) tabular export of the S3 block
for F, seed in ((0.01, -0.046 - 0.001j), (0.02, -0.0508977 - 0.0067223j)):
    picks = []
    for th in (0.18, 0.28):
        picks.append(nearest(F, 0.8, th, seed))
    e = picks[1]
    sh[F] = e.real - eyuk
    S3_ROWS.append({"case": f"ground F={F:.2f} lam=0.80", "F": F,
                    "lambda": 0.8, "Er": e.real, "Gamma": -2 * e.imag,
                    "shift": sh[F]})
    print(f"   F={F:.2f}: E = {e.real:+.7f}{e.imag:+.7f}j"
          f"   shift = {sh[F]*1e3:+.3f}e-3   (th-drift {abs(picks[1]-picks[0]):.0e})")
print(f"   ratio shift(0.02)/shift(0.01) = {sh[0.02]/sh[0.01]:.2f}"
      f"   (4.00 = quadratic)")

# (D-S4) export the field-halving datum.  `shift` = Er - E_Yukawa(lam=0.8)
# with the SAME F=0 anchor for both rows, so E_Yukawa = Er - shift is
# recoverable from either row; the manuscript's factor 2.24 is
# shift(F=0.02)/shift(F=0.01).  Values reported at theta_c = 0.28.
with (Path(__file__).with_name("d4_field_halving.csv")).open(
        "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["case", "F", "lambda", "Er",
                                       "Gamma", "shift"])
    w.writeheader()
    for r in S3_ROWS:
        w.writerow(r)

# ---------------- S4: tunneling exponent in its domain ------------------
print("\nS4 sub-barrier window at F=0.02 (F < Er^2/4 <=> |Er| > 0.283):")
D = Path(__file__).parent
M = {}
for row in csv.DictReader((D / "d3_map.csv").open()):
    M[(round(float(row["F"]), 3), round(float(row["lam"]), 3))] = (
        float(row["Er"]), float(row["Gamma"]))
for lam in (0.2, 0.25, 0.3):
    er = M[(0.02, lam)][0]
    print(f"   lam={lam}: Er={er:+.4f}  F_barrier=Er^2/4={er*er/4:.2e}"
          f"  {'sub-barrier' if 0.02 < er*er/4 else 'OVER-barrier'}")
ex = lambda er: np.exp(-(2 / (3 * 0.02)) * (2 * abs(er)) ** 1.5)
pred = ex(M[(0.02, 0.3)][0]) / ex(M[(0.02, 0.2)][0])
meas = M[(0.02, 0.3)][1] / M[(0.02, 0.2)][1]
print(f"   exponent ratio lam 0.2->0.3: predicted {pred:.0f}, measured {meas:.0f}")

print(f"\ntotal: {time.perf_counter()-t0:.1f} s")
