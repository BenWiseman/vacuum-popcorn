#!/usr/bin/env python3
"""Pair-production interaction length of UHE photons, computed from scratch.

Purpose: cross-check the numbers read off published figures in N_photon_horizon.md.
  1. CMB only (blackbody, T = 2.7255 K), Breit-Wheeler cross section.
  2. Radio background from the polynomial fit of Nitu et al. 2020 (arXiv:2004.13596,
     Eq. 39 and Table 2), alone and added to the CMB.
  3. Stress test with an ARCADE-2-like power law (A = 1.26 K at 1 GHz, index -2.6, values as
     quoted in arXiv:1110.5257), cut off below a chosen frequency. This is my own calculation,
     not a published curve.
Method: [lambda(E)]^-1 = 1/(8 E^2) Int d(eps) n(eps)/eps^2 Int_{4 me^2}^{4 eps E} s sigma(s) ds
(Protheroe and Biermann 1996, Eq. 23). Pure Python, no packages. Runs in under a second.

The CMB part is checked against four published numbers before anything is printed; the script
stops if any differs by more than 3 percent.
"""
import math, sys

me = 0.51099895e6            # eV
sigT = 6.6524587e-25         # cm^2
kB = 8.617333e-5             # eV/K
kT = kB * 2.7255
hbar_c = 1.973269804e-5      # eV cm
h_eVs = 4.135667696e-15      # eV s
c_cm = 2.99792458e10
Mpc = 3.0856775814913673e24  # cm
smin = 4 * me * me

def sigma_bw(s):
    x = 4 * me * me / s
    if x >= 1:
        return 0.0
    v = math.sqrt(1 - x)
    return 3.0 / 16.0 * sigT * (1 - v * v) * ((3 - v ** 4) * math.log((1 + v) / (1 - v)) - 2 * v * (2 - v * v))

def _G(smax):
    if smax <= smin:
        return 0.0
    n = 400
    a, b = math.log(smin), math.log(smax)
    h = (b - a) / n
    tot = 0.0
    for i in range(n + 1):
        s = math.exp(a + i * h)
        w = 1 if i in (0, n) else (4 if i % 2 else 2)
        tot += w * s * sigma_bw(s) * s
    return tot * h / 3

_N = 600
_lo, _hi = math.log(smin * 1.0000001), math.log(smin * 1e12)
_gx = [_lo + (_hi - _lo) * i / _N for i in range(_N + 1)]
_gy = [_G(math.exp(x)) for x in _gx]

def G(smax):
    l = math.log(smax)
    if l <= _gx[0]:
        return 0.0
    if l >= _gx[-1]:
        return _gy[-1]
    f = (l - _gx[0]) / (_gx[1] - _gx[0])
    i = int(f)
    r = f - i
    return _gy[i] * (1 - r) + _gy[i + 1] * r

def n_cmb(eps):
    x = eps / kT
    if x > 700:
        return 0.0
    return (eps ** 2 / (math.pi ** 2 * hbar_c ** 3)) / math.expm1(x)       # per eV per cm^3

_p = [-1.9847e+1, -2.9857e-1, -2.6984e-1, 9.5393e-2, -4.9059e-2, 4.4297e-3,
      7.6038e-3, -1.9690e-3, -2.2573e-4, 1.1762e-4, -9.9443e-6]            # Nitu et al. Table 2

def n_crb_nitu(eps):
    nu = eps / h_eVs
    x = math.log10(nu / 1e6)
    I = 10 ** sum(pi * x ** i for i, pi in enumerate(_p)) * 1e-4 * 6.241509074e18   # eV cm^-2 sr^-1 (per Hz s)
    return 4 * math.pi * I / (h_eVs * c_cm * eps)

def make_n_powerlaw(A, beta, nu_cut):
    kB_J, c_m = 1.380649e-23, 2.99792458e8
    def f(eps):
        nu = eps / h_eVs
        if nu < nu_cut:
            return 0.0
        T = A * (nu / 1e9) ** beta
        I = 2 * kB_J * T * nu ** 2 / c_m ** 2 * 1e-4 * 6.241509074e18
        return 4 * math.pi * I / (h_eVs * c_cm * eps)
    return f

def inv_lambda(E, nfun, lo, hi):
    a, b = math.log(max(smin / (4 * E), lo)), math.log(hi)
    n = 3000
    h = (b - a) / n
    tot = 0.0
    for i in range(n + 1):
        eps = math.exp(a + i * h)
        w = 1 if i in (0, n) else (4 if i % 2 else 2)
        tot += w * nfun(eps) / eps ** 2 * G(4 * eps * E) * eps
    return tot * h / 3 / (8 * E * E)                                          # per cm

def lam_cmb(E):
    return 1 / inv_lambda(E, n_cmb, 1e-9, 60 * kT) / Mpc

# ---- 1. validate the CMB part against published numbers (Mpc) ----
checks = [
    (2.2e15, 7.07e-3, "Vernetto and Lipari 2016 (arXiv:1608.01587): 7.07 kpc at 2.2 PeV"),
    (3.0e14, 72e-3,   "Vernetto and Lipari 2016: 72 kpc at 300 TeV"),
    (10 ** 17.3, 90e-3, "Auger 2017 (arXiv:1612.04155): 90 kpc at 10^17.3 eV"),
    (10 ** 18.5, 900e-3, "Auger 2017: 900 kpc at 10^18.5 eV"),
]
bad = False
for E, pub, label in checks:
    mine = lam_cmb(E)
    rel = abs(mine - pub) / pub
    flag = "ok" if rel < 0.03 else "FAIL"
    if rel >= 0.03:
        bad = True
    print("%-4s %-70s mine %.4g kpc  published %.4g kpc  (%.1f%%)" % (flag, label, mine * 1e3, pub * 1e3, rel * 100))
if bad:
    sys.exit("CMB check failed; do not use the numbers below")

# ---- 2. the table ----
print("\nlambda in Mpc (pair production, primary photon e-fold distance)")
print("%-9s %-10s %-12s %-14s %-22s" % ("E (eV)", "CMB only", "Nitu CRB", "CMB+Nitu CRB", "CMB+ARCADE-like>=1MHz"))
fl_nitu = (h_eVs * 1e3, h_eVs * 1e11)
arc = make_n_powerlaw(1.26, -2.6, 1e6)
for E in [1e19, 5e19, 1e20]:
    ic = inv_lambda(E, n_crb_nitu, *fl_nitu)
    im = inv_lambda(E, n_cmb, 1e-9, 60 * kT)
    ia = inv_lambda(E, arc, h_eVs * 1e6, h_eVs * 1e11)
    print("%-9.0e %-10.3g %-12.3g %-14.3g %-22.3g" % (E, 1 / im / Mpc, 1 / ic / Mpc, 1 / (ic + im) / Mpc, 1 / (ia + im) / Mpc))

# ---- 3. minimum of the CMB-only length ----
best = min(((lam_cmb(10 ** (14.5 + 0.05 * i)), 14.5 + 0.05 * i) for i in range(0, 25)))
print("\nCMB-only minimum: %.2f kpc at E = 10^%.2f eV" % (best[0] * 1e3, best[1]))

# ---- 4. loss over a 100 kpc path ----
for lam in (1.6, 2.4):
    print("lambda = %.1f Mpc: fraction absorbed over 100 kpc = %.1f%%, over 8.5 kpc = %.2f%%" %
          (lam, 100 * (1 - math.exp(-0.1 / lam)), 100 * (1 - math.exp(-0.0085 / lam))))
