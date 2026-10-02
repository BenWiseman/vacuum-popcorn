# p1_reference.py: the same wall encounters done naively at 100 significant digits (stdlib decimal),
# to check pk_lib.R's double-precision formulas. The particle is defined by (m_in, pt2, pn); its
# energy is recomputed here on the mass shell, as the R code assumes.
import csv, sys
from decimal import Decimal as Dm, getcontext
getcontext().prec = 100
path = sys.argv[1]
rows = list(csv.DictReader(open(path)))
worst = {"E": 0.0, "pn": 0.0, "A": 0.0}; mism = 0; border = 0; n = 0
for r in rows:
    m_in, m_out, pt2, pn, g = (Dm(float(r[k])) for k in ("m_in", "m_out", "pt2", "pn", "g"))
    E = (m_in*m_in + pt2 + pn*pn).sqrt()
    b = (1 - 1/(g*g)).sqrt()
    Ew = g*(E + b*pn); u = g*(pn + b*E)
    D2 = m_out*m_out - m_in*m_in
    if u <= 0: t, Ea, pa = 1, E, pn
    elif u*u < D2: t, Ea, pa = 2, g*(Ew + b*u), -g*(u + b*Ew)
    else:
        q = (u*u - D2).sqrt(); t, Ea, pa = 3, g*(Ew - b*q), g*(q - b*Ew)
    Aa = Ea + pa
    n += 1
    if t != int(float(r["type"])):
        # a disagreement only counts if it is not a borderline case (u^2 within 1e-12 of D2)
        if abs(u*u - D2) > Dm("1e-12")*D2: mism += 1
        else: border += 1
        continue
    if t == 1: continue
    for key, ref, got in (("E", Ea, r["E_after"]), ("pn", pa, r["pn_after"]), ("A", Aa, r["A_after"])):
        got = Dm(float(got))
        if ref != 0:
            worst[key] = max(worst[key], float(abs(got - ref)/abs(ref)))
print(f"cases {n}; type mismatches {mism}; borderline {border}")
print("largest relative error, double vs 100 digits: E_after %.2e, pn_after %.2e, light-cone A_after %.2e"
      % (worst["E"], worst["pn"], worst["A"]))
ok = mism == 0 and max(worst.values()) < 1e-9
print("PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
