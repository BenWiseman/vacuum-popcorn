# pk_lib.R: exact kinematics of a particle meeting a moving wall between two vacua.
#
# The particle has mass m_in on the side it comes from and m_out on the other side. The wall moves
# with Lorentz factor g along -n, where n is the unit normal pointing from the particle's side to the
# far side (for a collapsing pocket: n points outward and the wall moves inward). The particle's
# momentum splits into pn = p.n and a transverse part with square pt2, which the wall conserves.
# In the wall frame the energy is conserved; the particle is reflected while its normal momentum u
# satisfies u^2 < m_out^2 - m_in^2, and otherwise crosses with normal momentum sqrt(u^2 - D2).
# Every quantity is written so that nothing cancels at huge g or huge energy:
#   1 - beta = 1/(g^2 (1 + beta));  A = E + pn (the light-cone component the wall conserves as g -> oo).
encounter <- function(E, pn, pt2, m_in, m_out, g, A = NULL) {
  b <- sqrt(1 - 1/g^2); omb <- 1/(g^2*(1 + b))
  mT2 <- m_in^2 + pt2
  if (is.null(A)) A <- if (pn >= 0) E + pn else mT2/(E - pn)
  Ew <- g*(A - omb*pn)                        # wall-frame energy, g (E + b pn)
  u  <- g*(A - omb*E)                         # wall-frame normal momentum, g (pn + b E)
  if (u <= 0) return(list(type = "miss", E = E, pn = pn, A = A))
  D2 <- m_out^2 - m_in^2
  if (u^2 < D2) {                             # reflected: u -> -u
    Er <- g*(Ew + b*u); pr <- -g*(u + b*Ew)
    Ar <- mT2/((Ew + u)*g*(1 + b))            # E + pn after, without cancellation
    return(list(type = "reflect", E = Er, pn = pr, A = Ar, u = u))
  }
  q <- sqrt(u^2 - D2); MT2 <- m_out^2 + pt2   # crossed: now mass m_out
  Et <- g*(MT2/(Ew + q) + omb*q)
  pt <- g*(omb*Ew - MT2/(Ew + q))
  At <- (Ew + q)/(g*(1 + b))
  list(type = "transmit", E = Et, pn = pt, A = At, u = u)
}
# the limit g -> oo for a particle with light-cone component A and transverse momentum pt2
sweep_limit <- function(A, pt2, m_out) (A + (m_out^2 + pt2)/A)/2
