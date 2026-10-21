"""Plan D, Phase D3 v2: the (F, lambda) map with physical continuation.

v1 failed ("lost" at every lam>0) because the continuation window used a
fixed +0.09 guess while the true screening upshift is ~ +lam (the exact
first-order slope dE/dlam|_0 = <e^{-lam r}> = 1 for ANY normalized
state).  v2 predicts each point as

    E_pred(F, lam) = E_Yukawa(lam) + [E_prev - E_Yukawa(lam_prev)]

i.e. the exact F=0 bound-state curve (one cheap symmetric eigensolve
per lam, l=0 block carries the ground state) plus the slowly varying
complex field shift carried over from the previous point.  Steps are
fine (dlam=0.05), the search radius and the theta-agreement tolerance
scale with the running width, and a lost point retries once with a
larger rotation pair before being declared dissolved (tracking stops
for that F row after two consecutive losses).

D3b: n=2, m=0 is a DOUBLET (parabolic k=+-1), not a triplet -- v1's
label was wrong.  Both components are tracked separately in lam with
dlam=0.025; the blue (k=+1) component is followed through E_r=0 into
the positive-energy region.
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
print(f"[basis lmax=14, nb=70, DOF={15*70}, cert={S.cert:.0e}]")

THETAS = (0.18, 0.28)
THETAS_RETRY = (0.30, 0.42)


def track_point(F, lam, pred, radius, tol, thetas):
    """Eigenvalue nearest `pred` within `radius`, theta-stable to `tol`."""
    es = []
    for th in thetas:
        ev = np.linalg.eigvals(S.build(1.0, lam, F, th))
        ev = ev[(np.abs(ev - pred) < radius) & (ev.imag < 1e-8)]
        if not ev.size:
            return None
        es.append(ev[np.argmin(np.abs(ev - pred))])
    d = abs(es[1] - es[0])
    return (es[1], d) if d <= tol else None


# ---------------- Yukawa F=0 reference curve ----------------------------
LAMS = np.round(np.arange(0.0, 0.801, 0.05), 3)
print("\nYukawa F=0 ground reference (prediction anchor + physics column):")
eyuk = {}
for lam in LAMS:
    H = S.build(1.0, lam, 0.0, 0.0)
    eyuk[lam] = float(np.linalg.eigvalsh(0.5 * (H.real + H.real.T)).min())
print("   " + "  ".join(f"{l:.2f}:{eyuk[l]:+.4f}" for l in LAMS[::4]))

# ---------------- D3: ground-state map ----------------------------------
print("\nD3 v2 ground-state map:")
rows = []
for F in (0.02, 0.04, 0.06, 0.08, 0.10):
    prev, prev_lam, losses = None, None, 0
    for lam in LAMS:
        if losses >= 2:
            rows.append(dict(F=F, lam=lam, Er=np.nan, Gamma=np.nan,
                             Eyuk=eyuk[lam], qual=np.nan))
            continue
        if prev is None:
            pred, radius, tol = -0.52 + 0.0j, 0.08, 5e-3
        else:
            pred = eyuk[lam] + (prev - eyuk[prev_lam])
            radius = max(0.04, 4 * abs(prev.imag), 3 * abs(eyuk[lam] - eyuk[prev_lam]))
            tol = max(5e-3, 0.2 * abs(prev.imag))
        got = track_point(F, lam, pred, radius, tol, THETAS)
        if got is None:
            got = track_point(F, lam, pred, 1.5 * radius, 2 * tol, THETAS_RETRY)
        if got is None:
            losses += 1
            tag = "dissolved" if losses >= 2 else "lost"
            print(f"   F={F:.2f} lam={lam:.2f}: {tag}")
            rows.append(dict(F=F, lam=lam, Er=np.nan, Gamma=np.nan,
                             Eyuk=eyuk[lam], qual=np.nan))
            continue
        losses = 0
        e, q = got
        rows.append(dict(F=F, lam=lam, Er=e.real, Gamma=-2 * e.imag,
                         Eyuk=eyuk[lam], qual=q))
        print(f"   F={F:.2f} lam={lam:.2f}: E_r={e.real:+.7f}"
              f"  Gamma={-2*e.imag:.4e}  (q {q:.0e})")
        prev, prev_lam = e, lam
    print()

with (Path(__file__).with_name("d3_map.csv")).open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    for r in rows:
        w.writerow(r)

# ---------------- D3b: n=2 doublet vs lambda at F=0.02 ------------------
print("D3b n=2 doublet (m=0, parabolic k=-1 red / k=+1 blue) at F=0.02:")
man = []
comps = {"red": -0.2067 + 0.0j, "blue": -0.0890 + 0.0j}   # v1 lam=0 finds
state = {k: None for k in comps}
LAMS2 = np.round(np.arange(0.0, 0.301, 0.025), 3)
for lam in LAMS2:
    for name, seed in comps.items():
        if state[name] == "dead":
            continue
        prev = state[name]
        if prev is None:
            pred, radius, tol = seed, 0.05, 8e-3
        else:
            pred = prev[0] + (lam - prev[1]) * prev[2]
            radius = max(0.035, 3 * abs(prev[0].imag))
            tol = max(8e-3, 0.2 * abs(prev[0].imag))
        got = track_point(0.02, lam, pred, radius, tol, THETAS)
        if got is None:
            got = track_point(0.02, lam, pred, 1.5 * radius, 2 * tol,
                              THETAS_RETRY)
        if got is None:
            print(f"   lam={lam:.3f} {name}: lost -> stop")
            state[name] = "dead"
            continue
        e, q = got
        slope = 0.9 if prev is None else (e - prev[0]) / (lam - prev[1]) \
            if lam > prev[1] else 0.9
        state[name] = (e, lam, slope)
        man.append((name, lam, e.real, -2 * e.imag, q))
        print(f"   lam={lam:.3f} {name}: E_r={e.real:+.6f}"
              f"  Gamma={-2*e.imag:.3e}  (q {q:.0e})")

with (Path(__file__).with_name("d3_manifold.csv")).open("w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["comp", "lam", "Er", "Gamma", "qual"])
    for r in man:
        w.writerow(r)

print(f"\ntotal: {time.perf_counter()-t0:.1f} s")
