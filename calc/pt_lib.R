# pt_lib.R: follow one trapped particle, initially at rest, through a collapsing pocket.
#
# Thin wall, flat space, starting at rest at radius R0 = 1; vacuum energy eps above ours inside,
# tension sigma; k = eps R0 / (3 sigma). Energy conservation gives the wall Lorentz factor at radius x:
#   gamma(x) = (k (1 - x^3) + 1) / x^2,   so gamma - 1 = (1 - x)(k (1 + x + x^2) + 1 + x) / x^2.
# Units: c = 1, R0 = 1, energies in units of m_in (so m_in = 1 and m_out = mu).
# A particle at rest starts at radius x0. Every event happens at the wall. Between events the particle
# moves in a straight line through the centre (it only ever has radial momentum), so the geometry is
# one-dimensional and every timing question reduces to integrals of the wall's lag behind light,
# (1/beta - 1), and the particle's, (1/v - 1), computed as small numbers so nothing cancels.
source("pk_lib.R")
gam  <- function(x, k) (k*(1 - x^3) + 1)/x^2
gm1  <- function(x, k, omx = 1 - x) omx*(k*(1 + x + x^2) + 1 + x)/x^2
beta_x <- function(x, k, omx = 1 - x) { h <- gm1(x, k, omx); sqrt(h*(h + 2))/(h + 1) }
lag  <- function(x, k, omx = 1 - x) {          # 1/beta - 1 = 1/(g^2 (1 + b) b)
  g <- gm1(x, k, omx) + 1; b <- beta_x(x, k, omx); 1/(g^2*(1 + b)*b) }
I_lag <- function(a, b, k) {                    # integral of (1/beta - 1) from a to b, 0 <= a < b <= 1
  if (b <= a) return(0)
  tot <- 0; cut <- 0.5
  if (a < cut) { hi <- min(b, cut)
    tot <- tot + integrate(function(x) lag(x, k), a, hi, rel.tol = 1e-12, abs.tol = 0, subdivisions = 2000L)$value }
  if (b > cut) { lo <- max(a, cut)              # x = 1 - s^2 removes the 1/sqrt(1 - x) end point
    f <- function(s) lag(1 - s^2, k, s^2)*2*s
    tot <- tot + integrate(f, sqrt(1 - b), sqrt(1 - lo), rel.tol = 1e-12, abs.tol = 0, subdivisions = 2000L)$value }
  tot
}
# Every encounter obeys E_after < 2 g m_out (the wall-frame energy is below m_out for a reflection,
# and a crossing gives at most g m_out), checked at each step. Catch-ups from behind are counted:
# they are what lets a particle beat the head-on ceiling (2 m_out^2 - m_in^2)/m_in.
track <- function(x0, mu, k, xmin = 1e-6, max_events = 400) {
  E <- 1; p <- 0; A <- 1; x <- x0; dir <- "rest"; n <- 0; first <- NA; nc <- 0; gl <- NA; E1 <- NA
  res <- function(status, E) list(status = status, E = E, n = n, x = x, first = first, ncatch = nc, glast = gl, E1 = E1)
  repeat {
    n <- n + 1
    if (n > max_events) return(res("too_many", E))
    g <- gam(x, k); gl <- g
    pn <- switch(dir, rest = 0, out = p, inward = -p)
    if (dir == "inward") nc <- nc + 1
    r <- encounter(E, pn, 0, 1, mu, g, A = A)
    if (is.na(first)) first <- r$type
    if (r$type != "miss" && r$E > 2*g*mu*(1 + 1e-9)) stop("encounter bound E_after < 2 g m_out violated")
    if (r$type == "miss") return(res("miss", E))
    if (r$type == "transmit") return(res("escaped", r$E))
    # reflected: the particle now moves inward from radius x with energy E, light-cone A = E - p
    E <- r$E; p <- -r$pn; A <- r$A; delta <- A/p   # 1/v - 1
    if (is.na(E1)) E1 <- E                          # energy after the first reflection
    F0 <- I_lag(0, x, k) - x*delta                  # >0: particle reaches the centre first
    if (F0 > 0) {                                   # crosses the centre, meets the far side head-on
      G <- function(rho) F0 - 2*rho - I_lag(0, rho, k) - rho*delta
      rho <- uniroot(G, c(0, min(F0/2, x)), tol = 1e-14*F0)$root
      if (rho < xmin) { x <- rho; return(res("trapped_at_end", E)) }
      x <- rho; dir <- "out"; A <- E + p
    } else {                                        # the wall behind catches it up
      xpk <- uniroot(function(X) lag(X, k) - delta, c(1e-300^(1/8), x), tol = 1e-15)$root
      Ff <- function(X) I_lag(X, x, k) - (x - X)*delta
      xc <- uniroot(Ff, c(0, xpk), tol = 1e-15)$root
      if (xc < xmin) { x <- xc; return(res("trapped_at_end", E)) }
      x <- xc; dir <- "inward"
    }
  }
}
