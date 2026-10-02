#!/usr/bin/env python3
"""Scout M: check kinematic facts (a), (b), (c) by explicit Lorentz boosts, c = 1.

Setup: plasma frame, wall moves along +z with Lorentz factor gamma. The particle approaches from
the plasma side, so in the wall frame its normal momentum is -u (u > 0). The wall frame conserves
energy and p_perp. Mass is m_in before the wall and m_out after it. Transmitted if
E_w^2 - p_perp^2 - m_out^2 >= 0, otherwise reflected (normal momentum flips to +u).

Checks, all in 60-digit Decimal arithmetic:
 (a) particle at rest, gamma -> large: E_after -> (m_out^2 + m_in^2)/(2 m_in), and E - p_z is kept.
 (b) particle at rest, reflected: E_after = gamma^2 (1 + beta^2) m_in.
 (c1) every reflection: E_before * E_after = gamma^2 (m_in^2 + p_perp^2) + u^2.
      The form with p_perp dropped must FAIL for p_perp != 0 (this shows the check can fail).
 (c2) one bounce from rest: E_after tends to (2 m_out^2 - m_in^2)/m_in as gamma -> m_out/m_in from below;
      just above that gamma the particle is transmitted.
Run: python3 M_kin_check.py   (prints PASS/FAIL lines, exits 1 on any failure)
"""
import random
import sys
from decimal import Decimal as D, getcontext

getcontext().prec = 60
FAIL = []


def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + (("  " + detail) if detail else ""))
    if not ok:
        FAIL.append(name)


def boost(E, pz, beta):
    g = 1 / (1 - beta**2).sqrt()
    return g * (E - beta * pz), g * (pz - beta * E)


def unboost(E, pz, beta):
    g = 1 / (1 - beta**2).sqrt()
    return g * (E + beta * pz), g * (pz + beta * E)


def cross(E, pz, ptsq, m_out, beta):
    """Return ('R'|'T', E_after, pz_after, u) in the plasma frame."""
    Ew, pzw = boost(E, pz, beta)
    assert pzw < 0, "particle must approach the wall"
    u = -pzw
    pz2sq = Ew**2 - ptsq - m_out**2
    if pz2sq < 0:
        Ea, pza = unboost(Ew, +u, beta)
        return "R", Ea, pza, u
    Ea, pza = unboost(Ew, -pz2sq.sqrt(), beta)
    return "T", Ea, pza, u


# (b) at rest, reflected
m_in, m_out = D(1), D(10)
for g in (D(2), D(5), D("9.9")):
    beta = (1 - 1 / g**2).sqrt()
    st, Ea, pza, u = cross(m_in, D(0), D(0), m_out, beta)
    check("b: at-rest reflection gamma=%s" % g, st == "R" and abs(Ea - g**2 * (1 + beta**2) * m_in) < D("1e-40"),
          "E_after/m_in=%.6f" % float(Ea / m_in))

# (a) at rest, transmitted, gamma large
for (mi, mo) in ((D(1), D(10)), (D("0.5"), D(3))):
    for g in (D(10) ** 6, D(10) ** 9):
        beta = (1 - 1 / g**2).sqrt()
        st, Ea, pza, u = cross(mi, D(0), D(0), mo, beta)
        target = (mo**2 + mi**2) / (2 * mi)
        check("a: m_in=%s m_out=%s gamma=%.0e" % (mi, mo, g), st == "T" and abs(Ea - target) / target < D("1e-10")
              and abs((Ea - pza) - mi) / mi < D("1e-10"), "E_after=%.10f target=%.10f" % (float(Ea), float(target)))

# (c1) product relation, random reflections
random.seed(1)
dev_T = D(0)
dev_0 = D(0)
n = 0
for _ in range(20000):
    g = D(random.uniform(1.05, 50))
    beta = (1 - 1 / g**2).sqrt()
    mi = D(random.uniform(0.1, 3))
    mo = mi + D(random.uniform(0.1, 30))
    pt = D(random.uniform(0.5, 5))
    pz = D(random.uniform(-40, 40))
    E = (mi**2 + pt**2 + pz**2).sqrt()
    if boost(E, pz, beta)[1] >= 0:
        continue
    st, Ea, pza, u = cross(E, pz, pt**2, mo, beta)
    if st != "R":
        continue
    n += 1
    dev_T = max(dev_T, abs(E * Ea - (g**2 * (mi**2 + pt**2) + u**2)) / (E * Ea))
    dev_0 = max(dev_0, abs(E * Ea - (g**2 * mi**2 + u**2)) / (E * Ea))
check("c1: E*E' = gamma^2 (m_in^2 + pt^2) + u^2 (n=%d reflections, pt>0)" % n, dev_T < D("1e-40"), "max rel dev %.2e" % float(dev_T))
check("c1 control: form without pt FAILS when pt>0 (check is sensitive)", dev_0 > D("0.1"), "max rel dev %.3f" % float(dev_0))

# (c1) pt = 0 only: the form in the letter holds exactly
dev = D(0)
n0 = 0
for _ in range(20000):
    g = D(random.uniform(1.05, 50))
    beta = (1 - 1 / g**2).sqrt()
    mi = D(random.uniform(0.1, 3))
    mo = mi + D(random.uniform(0.1, 30))
    pz = D(random.uniform(-40, 40))
    E = (mi**2 + pz**2).sqrt()
    if boost(E, pz, beta)[1] >= 0:
        continue
    st, Ea, pza, u = cross(E, pz, D(0), mo, beta)
    if st != "R":
        continue
    n0 += 1
    dev = max(dev, abs(E * Ea - (g**2 * mi**2 + u**2)) / (E * Ea))
check("c1: pt=0, E*E' = gamma^2 m_in^2 + u^2 (n=%d)" % n0, dev < D("1e-40"), "max rel dev %.2e" % float(dev))

# (c2) one-bounce ceiling from rest
for (mi, mo) in ((D(1), D(10)), (D("0.5"), D(3)), (D(1), D("1.5"))):
    gc = mo / mi
    ceil = (2 * mo**2 - mi**2) / mi
    g = gc * (1 - D("1e-12"))
    beta = (1 - 1 / g**2).sqrt()
    st, Ea, _, _ = cross(mi, D(0), D(0), mo, beta)
    check("c2: below gamma_c, m_in=%s m_out=%s reflects, E_after -> %s" % (mi, mo, ceil),
          st == "R" and abs(Ea - ceil) / ceil < D("1e-9"), "E_after=%.8f" % float(Ea))
    g = gc * (1 + D("1e-6"))
    beta = (1 - 1 / g**2).sqrt()
    st, Ea, _, _ = cross(mi, D(0), D(0), mo, beta)
    check("c2: above gamma_c transmits", st == "T")

print("ALL PASS" if not FAIL else "FAILED: %s" % FAIL)
sys.exit(1 if FAIL else 0)
