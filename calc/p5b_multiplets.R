# p5b_multiplets.R: at the photon cap of p5, how often would a pocket population give one event, and
# how often two or more at once, in a 3000 km^2 array? Pockets of mass Mp make up all of the Galaxy's
# dark matter (NFW as in p5); each pop puts the fraction f = 1 of its mass-energy into chi at
# E_star = 1e11 GeV; each chi gives one detectable particle (an upper estimate); events per pop are
# Poisson with mean m(d) = N_chi A / (4 pi d^2).
kpc <- 3.0857e21; yr <- 3.15576e7; GeV_g <- 1.78266e-24; E_star <- 1e11; A <- 3000*1e10
rs <- 20; Rsun <- 8.2; nfw <- function(r) 1/((r/rs)*(1 + r/rs)^2); rho_s <- 0.4/nfw(Rsun)  # GeV/cm^3
fG <- 1.6e-32                                    # per s, the photon cap from p5
cat("pocket mass | pops/yr with >=1 event | with >=2 | with >=3 | events/yr in all | share of events in groups\n")
for (Mg in c(1e15, 1e18, 1e20, 1e22, 1e24)) {
  Mp <- Mg/GeV_g; Nchi <- Mp/E_star                # pocket mass in GeV, chi per pop
  # integrate over distance d from the Sun and direction (cos angle from the Galactic centre)
  cb <- seq(-1, 1, length.out = 401); cb <- (cb[-1] + cb[-length(cb)])/2     # midpoints, never through r = 0
  shell <- function(d) sapply(d, function(dd) mean(rho_s*nfw(sqrt(Rsun^2 + dd^2 - 2*Rsun*dd*cb)))*4*pi*(dd*kpc)^2*kpc)
  # shell(d): GeV of dark matter per kpc of distance at distance d (integral of rho over the shell)
  f_d <- function(d, n) {                          # pops per year per kpc of distance, with >= n events
    m <- Nchi*A/(4*pi*(d*kpc)^2)
    p <- if (n == 1) -expm1(-m) else ppois(n - 1, m, lower.tail = FALSE)
    shell(d)/Mp*fG*yr*p
  }
  segs <- c(1e-9, 1e-3, 1, 7, 8.2, 9.5, 30, 300)
  integ <- function(fun, ...) sum(sapply(seq_len(length(segs) - 1), function(i)
    integrate(fun, segs[i], segs[i + 1], ..., subdivisions = 4000L, rel.tol = 1e-7)$value))
  R <- sapply(1:3, function(n) integ(f_d, n = n))
  ev <- integ(function(d) shell(d)/Mp*fG*yr*Nchi*A/(4*pi*(d*kpc)^2))
  grp <- 1 - R[1]*1/ev                              # share of events not arriving as singles (approx.)
  cat(sprintf("  %.0e g   | %10.3g            | %8.3g | %8.3g | %10.3g      | %.3f\n", Mg, R[1], R[2], R[3], ev,
              max(0, 1 - (R[1] - R[2])/ev)))
}
cat("\n(events/yr in all is independent of pocket mass at fixed f Gamma: the photon cap fixes the mean flux)\n")
