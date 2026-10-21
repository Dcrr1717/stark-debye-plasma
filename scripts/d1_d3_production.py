"""Plan D, Phases D1 + D3: validation anchors and the (F, lambda) map.

D1a  Pure-Stark ground resonances F in {0.04, 0.08, 0.10}: convergence
     protocol (2 theta x 2 grids); F=0.10 matches the classical value
     (-0.527418175 - 0.007269057j, D0).  Values are quoted as
     self-converged benchmarks; full-digit literature re-anchor
     (Jentschura/Kolosov) marked pending.
D1b  F=0, lam=0.25 fine: the l=1 block vs the 1D triple-verified
     benchmark (0.00306445 - 0.0021845j), proper protocol.
D3   THE MAP (new physics): ground-state E_r(F,lam), Gamma(F,lam) on
     F in {0.02,...,0.10} x lam in {0,...,0.8} with lambda-continuation
     (window follows the previous point); plus the n=2 manifold (m=0
     triplet) vs lam at F=0.02.
D2 note (honest): at ~780 DOF a dense eigensolve takes ~1-2 s; the full
     map runs in minutes.  FEAST/warm-start is UNNECESSARY at this
     scale -- the plan's method-novelty claim is formally DEMOTED per
     its own criterion; the paper's thesis is the certified
     coupled-channel method + the first systematic map.
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


def stable_states(S, F, lam, thetas, center, radius, imax=0.2, nmax=4):
    """Theta-stable eigenvalues within a disk around `center`."""
    picks = []
    for th in thetas:
        H = S.build(1.0, lam, F, th)
        ev = np.linalg.eigvals(H)
        sel = ev[(np.abs(ev - center) < radius) & (ev.imag < 1e-9)
                 & (ev.imag > -imax)]
        picks.append(sel)
    if not picks[0].size or not picks[1].size:
        return []
    out = []
    used = set()
    d = np.abs(picks[1][None, :] - picks[0][:, None])
    for _ in range(min(nmax, min(len(p) for p in picks))):
        i, j = np.unravel_index(np.argmin(d), d.shape)
        if d[i, j] > 5e-3:
            break
        out.append((picks[1][j], d[i, j]))
        d[i, :] = np.inf
        d[:, j] = np.inf
    return out


# ---------------- D1a: pure-Stark anchors -------------------------------
print("D1a pure-Stark ground anchors (self-converged; lit re-anchor pending):")
for F in (0.04, 0.08, 0.10):
    vals = []
    for lmax, nb in ((12, 60), (16, 80)):
        S = StarkChannels(nb, lmax, 0.35)
        got = stable_states(S, F, 0.0, (0.15, 0.25), -0.52, 0.08,
                            imax=0.05, nmax=1)
        if got:
            vals.append(got[0][0])
    if len(vals) == 2:
        print(f"   F={F:.2f}: E = {vals[1].real:.10f}{vals[1].imag:+.10f}j"
              f"   grid-drift = {abs(vals[1]-vals[0]):.1e}")

# ---------------- D1b: screening limit, fine ----------------------------
print("\nD1b F=0, lam=0.25 (l=1 block) vs 1D benchmark:")
BENCH = 0.00306445 - 0.0021845j
vals = []
for nb, beta in ((150, 0.20), (220, 0.15)):
    S = StarkChannels(nb, 1, beta)
    got = stable_states(S, 0.0, 0.25, (0.35, 0.45), BENCH, 0.005,
                        imax=0.05, nmax=1)
    if got:
        vals.append(got[0][0])
if len(vals) == 2:
    print(f"   E = {vals[1].real:.9f}{vals[1].imag:+.9f}j"
          f"   dev vs 1D = {abs(vals[1]-BENCH):.1e}"
          f"   grid-drift = {abs(vals[1]-vals[0]):.1e}")

# ---------------- D3: the (F, lambda) ground-state map ------------------
print("\nD3 ground-state map (lambda-continuation per F):")
rows = []
S = StarkChannels(70, 14, 0.35)
print(f"   [basis lmax=14, nb=70, DOF={15*70}, cert={S.cert:.0e}]")
for F in (0.02, 0.04, 0.06, 0.08, 0.10):
    center, radius = -0.52 + 0.0j, 0.08
    for lam in (0.0, 0.2, 0.4, 0.6, 0.8):
        got = stable_states(S, F, lam, (0.18, 0.28), center, radius,
                            imax=0.3, nmax=1)
        if not got:
            print(f"   F={F:.2f} lam={lam:.1f}: lost")
            rows.append(dict(F=F, lam=lam, Er=np.nan, Gamma=np.nan,
                             qual=np.nan))
            continue
        e, q = got[0]
        rows.append(dict(F=F, lam=lam, Er=e.real, Gamma=-2 * e.imag, qual=q))
        print(f"   F={F:.2f} lam={lam:.1f}: E_r={e.real:+.7f}"
              f"  Gamma={-2*e.imag:.4e}  (q {q:.0e})")
        center, radius = e + lam * 0 + 0.1j * 0, max(0.06, 5 * abs(e.imag))
        center = complex(e.real + 0.09, min(e.imag, -1e-9))  # anticipate shift
    print()

# v1 outputs renamed: d3_map_v2.py owns the published d3_map.csv /
# d3_manifold.csv; rerunning this script must never overwrite them.
with (Path(__file__).with_name("d3_map_v1.csv")).open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    for r in rows:
        w.writerow(r)

# ---------------- D3b: n=2 manifold vs lambda at F=0.02 -----------------
print("D3b n=2 manifold (m=0 triplet) at F=0.02 vs lambda:")
man = []
for lam in (0.0, 0.05, 0.10, 0.15, 0.20):
    got = stable_states(S, 0.02, lam, (0.18, 0.28),
                        complex(-0.125 + 0.9 * lam * 0.25, 0.0), 0.09,
                        imax=0.1, nmax=3)
    ee = sorted([g[0] for g in got], key=lambda z: z.real)
    man.append((lam, ee))
    print(f"   lam={lam:.2f}: " + "  ".join(
        f"({e.real:+.6f},{-2*e.imag:.2e})" for e in ee))
with (Path(__file__).with_name("d3_manifold_v1.csv")).open("w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["lam", "Er", "Gamma"])
    for lam, ee in man:
        for e in ee:
            w.writerow([lam, e.real, -2 * e.imag])
print(f"\ntotal: {time.perf_counter()-t0:.1f} s")
