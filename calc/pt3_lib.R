# pt3_lib.R: follow one trapped particle with any initial momentum through a collapsing pocket, in 3D.
# Same wall and units as pt_lib.R (R0 = 1, c = 1, m_in = 1, m_out = mu). Between events the particle
# moves on a straight chord. The wall reaches radius x at time T(x) = (1 - x) + Ilag(x, 1), with
# Ilag(a, b) the integral of (1/beta - 1) from a to b, tabulated once per k so that every timing
# question is a difference of small, accurately known numbers.
source("pt_lib.R")
make_lagtab <- function(k) {
  # lower part: C(x) = Ilag(0, x) on a log grid in x (C ~ x^5 for small x); upper part: L(s) =
  # Ilag(1 - s^2, 1) on a grid in s = sqrt(1 - x), which removes the 1/sqrt(1 - x) end point
  xl <- 10^seq(-40, log10(0.5), length.out = 3000)
  seg <- vapply(seq_along(xl), function(i) if (i == 1) 0 else
    integrate(function(x) lag(x, k), xl[i - 1], xl[i], rel.tol = 1e-13, abs.tol = 0)$value, 0)
  C0 <- (k + 1)^(-2)*xl[1]^5/10                  # integral from 0 to the first grid point, x^4/(2(k+1)^2)
  Cl <- C0 + cumsum(seg)
  su <- seq(0, sqrt(0.5), length.out = 3000)
  f <- function(s) lag(1 - s^2, k, s^2)*2*s
  segu <- vapply(seq_along(su), function(i) if (i == 1) 0 else
    integrate(f, su[i - 1], su[i], rel.tol = 1e-13, abs.tol = 0)$value, 0)
  Lu <- cumsum(segu)                             # Ilag(1 - s^2, 1)
  lc <- splinefun(log(xl), log(Cl), method = "monoH.FC")
  lu <- splinefun(su, Lu, method = "monoH.FC")
  Chalf <- Cl[length(Cl)]; Lhalf <- Lu[length(Lu)]
  Cfun <- function(x) {                          # Ilag(0, x) for 0 <= x <= 1
    out <- numeric(length(x))
    lo <- x <= 0.5 & x > xl[1]; out[lo] <- exp(lc(log(x[lo])))
    tiny <- x <= xl[1]; out[tiny] <- (k + 1)^(-2)*x[tiny]^5/10
    hi <- x > 0.5; out[hi] <- Chalf + Lhalf - lu(sqrt(1 - x[hi]))
    out }
  list(C = Cfun, k = k)
}
Ilag_t <- function(tab, a, b) tab$C(b) - tab$C(a)   # signed

next_hit <- function(r, v, E, m2, A_dep, Xref, at_wall, tab, xmin) {
  # first s > 0 with T(|r + v s|) = T(Xref) + s, i.e. G(s) = (Xref - rho) + Ilag(rho, Xref) - s = 0.
  # at_wall: the particle starts on the wall (|r| = Xref) after an event; A_dep = E + p.n there.
  rr <- sum(r*r); rv <- sum(r*v); v2 <- 1 - m2/E^2
  Dfun <- function(s) {                           # Xref - s - rho(s), without cancellation
    rho <- sqrt(max(rr + 2*s*rv + s*s*v2, 0))
    den <- Xref - s + rho
    if (den > 0.25*(Xref + rho)) {
      a1 <- if (at_wall) Xref*A_dep/E else Xref + rv     # Xref + r.v
      num <- (Xref^2 - rr)*(!at_wall) - 2*s*a1 + s*s*m2/E^2
      c(num/den, rho)
    } else c(Xref - s - rho, rho)
  }
  G <- function(s) { d <- Dfun(s); d[1] + Ilag_t(tab, min(d[2], 1), Xref) }
  # time left until the wall reaches xmin
  s_end <- (Xref - xmin) + Ilag_t(tab, xmin, Xref)
  if (G(s_end) > 0) return(list(hit = FALSE))
  # scan for the first sign change on a log grid
  s0 <- if (at_wall) 1e-14*Xref else 0
  grid <- c(s0, s0 + (s_end - s0)*10^seq(-15, 0, length.out = 90))
  gv <- vapply(grid, G, 0)
  j <- which(gv[-1] <= 0 & gv[-length(gv)] > 0)[1]
  if (is.na(j)) return(list(hit = FALSE))
  s <- uniroot(G, c(grid[j], grid[j + 1]), tol = 1e-15*max(grid[j + 1], 1e-300))$root
  list(hit = TRUE, s = s, rh = r + v*s)
}

track3 <- function(r0, p0, mu, tab, xmin = 1e-12, max_events = 400) {
  k <- tab$k; m2 <- 1; E <- sqrt(1 + sum(p0^2)); p <- p0; r <- r0
  Xref <- 1; at_wall <- FALSE; A_dep <- NA; n <- 0; nc <- 0
  repeat {
    n <- n + 1
    if (n > max_events) return(list(status = "too_many", E = E, n = n, ncatch = nc))
    v <- p/E
    h <- next_hit(r, v, E, m2, A_dep, Xref, at_wall, tab, xmin)
    if (!h$hit) return(list(status = "trapped_at_end", E = E, n = n, ncatch = nc))
    rh <- h$rh; x <- sqrt(sum(rh*rh)); nrm <- rh/x
    if (x < xmin) return(list(status = "trapped_at_end", E = E, n = n, ncatch = nc))
    g <- gam(x, k)
    pn <- sum(p*nrm); cr <- c(p[2]*nrm[3] - p[3]*nrm[2], p[3]*nrm[1] - p[1]*nrm[3], p[1]*nrm[2] - p[2]*nrm[1])
    pt2 <- sum(cr*cr); ptv <- p - pn*nrm
    A <- if (pn >= 0) E + pn else (m2 + pt2)/(E - pn)
    if (pn < 0) nc <- nc + 1
    e <- encounter(E, pn, pt2, 1, mu, g, A = A)
    if (e$type != "miss" && e$E > 2*g*sqrt(mu^2 + pt2)*(1 + 1e-9)) stop("encounter bound violated")
    if (e$type == "transmit") return(list(status = "escaped", E = e$E, n = n, ncatch = nc))
    if (e$type == "miss") return(list(status = "miss", E = E, n = n, ncatch = nc))
    E <- e$E; p <- e$pn*nrm + ptv; r <- rh; Xref <- x; at_wall <- TRUE; A_dep <- e$A
  }
}
