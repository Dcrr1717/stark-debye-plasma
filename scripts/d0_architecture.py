"""Plan D, Gate D0: architecture benchmark -- coupled-channel exact
basis for the non-separable Stark problem.

The Stark term F z couples only l <-> l+1 (m = 0 sector), so the 2D
problem is a BLOCK-TRIDIAGONAL radial system: per channel our exact
Laguerre basis (shared beta and quadrature grid across channels so the
cross-l dipole elements <phi_i^l| r |phi_j^l'> are exact contractions),
complex-scaled analytically:

  H(theta) = diag_l [ e^{-2i th} T_l + V_th,l ]
           + F e^{+i th} sum_l c_l ( |l><l+1| x R + h.c. ),
  c_l = (l+1)/sqrt((2l+1)(2l+3))          (m = 0 Gaunt factor).

m != 0 (added for editorial task D-N1): m is conserved by F z, so the
same block-tridiagonal structure holds in any m sector; only the
channel list (l = |m| .. lmax) and the Gaunt factor change,

  c_l(m) = sqrt((l+1-m)(l+1+m)) / sqrt((2l+1)(2l+3)),

which is the FACTORED form actually implemented in `cl()` and is what
reduces to the expression above bit-for-bit at m = 0.  The algebraically
equivalent sqrt(((l+1)^2 - m^2)/((2l+1)(2l+3))) is 1 ulp off at
l = 0, 3, 5, 6, 9, 13 and must not be substituted here.  `StarkChannels` takes
an OPTIONAL `m=0` keyword; the default path is untouched and produces
bit-identical matrices (see d6_m1_probe.py, block R1).

PRE-REGISTERED PREDICTIONS:
P-D0.1 Hydrogen ground Stark resonance at F=0.1 a.u.: theta-stable
       value converging geometrically in l_max and nb, near the
       classical complex-rotation value E ~ -0.527418 - 0.00727i
       (memory-level anchor; to be re-anchored against Jentschura/
       Kolosov tables in D1 -- the gate here is INTERNAL convergence).
P-D0.2 F=0, lam=0.25 block-diagonal limit: the l=1 block reproduces
       the 1D engine's triple-verified 2p Yukawa resonance
       (0.00306445 - 0.0021845i) -- machinery self-check.
P-D0.3 DOF efficiency: convergence with N = (l_max+1)*nb <~ 1000 DOF,
       i.e. orders of magnitude below the 2D-mesh requirement --
       architecture (ii) takes the gate; the mesh becomes the D1
       cross-check tool.
"""

import sys
import time
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from basis import _q_functions  # noqa: E402

t0 = time.perf_counter()


class StarkChannels:
    """Coupled-channel complex-scaled Stark + screened Coulomb."""

    def __init__(self, nb, lmax, beta, nquad=None, m=0):
        self.nb, self.lmax, self.beta = nb, lmax, beta
        self.m = abs(int(m))
        # channels present in this m sector: l = |m| .. lmax
        # (at m = 0 this is range(lmax+1), i.e. the original enumeration)
        self.ls = list(range(self.m, lmax + 1))
        if not self.ls:
            raise ValueError(f"no channels: |m|={self.m} > lmax={lmax}")
        xmax = 4.0 * nb + 24.0 * (lmax + 1) + 120.0
        nq = nquad or max(3000, 16 * nb)
        t, w = leggauss(nq)
        self.x = 0.5 * xmax * (t + 1.0)
        self.w = 0.5 * xmax * w
        self.q, self.dq = {}, {}
        for l in self.ls:
            self.q[l], self.dq[l] = _q_functions(nb, l, self.x)
        # certification: worst overlap error across channels
        self.cert = max(
            float(np.max(np.abs((self.q[l] * self.w) @ self.q[l].T
                                - np.eye(nb))))
            for l in self.ls)

    def cl(self, l):
        """Dipole (Gaunt) factor <l+1,m| cos | l,m> for this m sector.

        sqrt(((l+1)^2 - m^2)/((2l+1)(2l+3))), written factored so that at
        m = 0 the numerator sqrt((l+1)*(l+1)) is EXACTLY l+1 in IEEE-754
        and the value is bit-identical to the original
        (l+1)/sqrt((2l+1)(2l+3)).  (The unfactored form
        sqrt(((l+1)**2 - m**2)/((2l+1)(2l+3))) is 1 ulp off at
        l = 0,3,5,6,9,13,... because sqrt(a^2/b) != a/sqrt(b) under
        round-to-nearest.)
        """
        m = self.m
        return (np.sqrt((l + 1 - m) * (l + 1 + m))
                / np.sqrt((2 * l + 1) * (2 * l + 3)))

    def build(self, Z, lam, F, theta):
        nb, lmax, b = self.nb, self.lmax, self.beta
        ls = self.ls
        N = len(ls) * nb
        H = np.zeros((N, N), dtype=complex)
        ph1, ph2 = np.exp(-1j * theta), np.exp(-2j * theta)
        r = self.x / (2 * b)
        for k, l in enumerate(ls):
            qw = self.q[l] * self.w
            T = 2 * b * b * ((self.dq[l] * self.w) @ self.dq[l].T)
            if l > 0:
                T += 2 * b * b * l * (l + 1) * ((qw / self.x**2) @ self.q[l].T)
            T = 0.5 * (T + T.T)                # symmetrize; T already
                                               # includes the 1/2 of p^2/2
            if lam == 0.0:
                V = -Z * 2 * b * ((qw / self.x) @ self.q[l].T)
                Vth = ph1 * 0.5 * (V + V.T)
            else:
                f = np.exp(-lam * r * np.exp(1j * theta))
                V = -Z * 2 * b * ((qw * f / self.x) @ self.q[l].T)
                Vth = ph1 * 0.5 * (V + V.T)
            i0 = k * nb
            H[i0:i0 + nb, i0:i0 + nb] = ph2 * T + Vth
        # Stark coupling F r cos(theta_ang), scaled r -> r e^{i theta}
        for k, l in enumerate(ls[:-1]):
            cl = self.cl(l)
            R = ((self.q[l] * self.w * r) @ self.q[l + 1].T)
            blk = F * cl * np.exp(1j * theta) * R
            i0, j0 = k * nb, (k + 1) * nb
            H[i0:i0 + nb, j0:j0 + nb] = blk
            H[j0:j0 + nb, i0:i0 + nb] = blk.T
        return H


def stark_resonance(nb, lmax, beta, F, thetas=(0.20, 0.30),
                    win=((-0.60, -0.40), (-0.10, -1e-7))):
    S = StarkChannels(nb, lmax, beta)
    picks = []
    for th in thetas:
        H = S.build(1.0, 0.0, F, th)
        ev = np.linalg.eigvals(H)
        sel = ev[(ev.real > win[0][0]) & (ev.real < win[0][1])
                 & (ev.imag > win[1][0]) & (ev.imag < win[1][1])]
        if sel.size:
            picks.append(sel)
    if len(picks) < 2:
        return None, np.inf, S.cert
    d = np.abs(picks[1][None, :] - picks[0][:, None])
    i, j = np.unravel_index(np.argmin(d), d.shape)
    return picks[1][j], d[i, j], S.cert


# ---------------- P-D0.0: hydrogen-levels certification -----------------
# (added after this exact check caught a halved-kinetic factor bug that
# the window-based P-D0.2 pick had coincidentally masked)
S0 = StarkChannels(60, 2, 0.35)
H0 = S0.build(1.0, 0.0, 0.0, 0.0)
ev0 = np.sort(np.linalg.eigvals(H0).real)
exact = [-0.5, -0.125, -0.125, -1.0 / 18, -1.0 / 18]
dev0 = max(abs(a - b) for a, b in zip(ev0[:5], exact))
print(f"P-D0.0 pure-hydrogen certification: lowest5 dev = {dev0:.1e}"
      + ("  OK" if dev0 < 1e-9 else "  *** FAIL ***"))
if not (dev0 < 1e-9):
    raise RuntimeError(
        f"P-D0.0 hydrogen certification FAILED (dev={dev0:.1e}) -- "
        "aborting: no production run may proceed on a failed assembly")

# ---------------- P-D0.2: block-diagonal self-check ---------------------
print("\nP-D0.2 F=0, lam=0.25 limit (l=1 block = 1D engine):")
S = StarkChannels(60, 2, 0.25)
for th in (0.35,):
    H = S.build(1.0, 0.25, 0.0, th)
    ev = np.linalg.eigvals(H)
    sel = ev[(ev.real > 1e-3) & (ev.real < 0.02)
             & (ev.imag < -1e-4) & (ev.imag > -0.05)]
    e = sel[np.argmin(np.abs(sel - (0.00306445 - 0.0021845j)))]
    print(f"   theta={th}: {e:.8f}   (1D benchmark 0.00306445-0.00218450j)"
          f"   dev={abs(e - (0.00306445 - 0.0021845j)):.1e}"
          f"   cert={S.cert:.0e}")

# ---------------- P-D0.1 / P-D0.3: Stark F=0.1 convergence --------------
print("\nP-D0.1 hydrogen Stark F=0.1: convergence in l_max and nb:")
print(f"{'lmax':>5} {'nb':>4} {'DOF':>5} {'E (theta-stable)':>26}"
      f" {'th-stab':>8} {'cert':>7}")
ref = None
for lmax, nb in ((4, 40), (8, 40), (12, 40), (12, 60), (16, 60)):
    e, q, cert = stark_resonance(nb, lmax, 0.35, 0.1)
    if e is None:
        print(f"{lmax:>5} {nb:>4} {(lmax+1)*nb:>5}   none")
        continue
    tag = "" if ref is None else f"   d_prev={abs(e - ref):.1e}"
    print(f"{lmax:>5} {nb:>4} {(lmax+1)*nb:>5} {e.real:13.9f}{e.imag:+13.9f}j"
          f" {q:>8.0e} {cert:>7.0e}{tag}")
    ref = e
print(f"\n   memory-level anchor: -0.527418 - 0.007269j (re-anchor in D1)")
print(f"total: {time.perf_counter()-t0:.1f} s")
