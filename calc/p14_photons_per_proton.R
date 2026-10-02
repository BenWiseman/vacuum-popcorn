# p14_photons_per_proton.R: how many photons come with each proton when a pop's escapee decays to quarks.
#
# The letter's caps on pop protons (p13) divide each photon limit by this ratio. The figure usually
# quoted, 2 to 3, is from Aloisio, Berezinsky and Kachelriess (hep-ph/0307279): "At x ~ 1e-3 this
# ratio is characterized by a value of 2 - 3 only", for a particle of 1e14 GeV decaying at rest, x being
# the share of its energy that a product carries. A pop's escapee differs in two ways. It is moving,
# with a Lorentz factor of 1e4 or more, and a proton seen near 1e20 eV carries a large share of its
# energy, not a thousandth.
#
# Their equations (18) and (19) give the spectra from the hadron spectrum D_h(x):
#   nucleons  D_N(x) = eps_N D_h(x),              eps_N  = 0.12 +- 0.02
#   photons   D_g(x) = (2/3) eps_pi I(x),         eps_pi = 0.73 +- 0.03,   I(x) = int_x^1 D_h(y) dy/y
# (a neutral pion is one pion in three and gives two photons, each uniform in energy up to the pion's).
# A product emitted in a random direction by a fast parent has a lab energy uniform between 0 and its
# share x of the parent's energy E_X, so in the lab every spectrum is smeared once more by int dy/y:
#   nucleons  N(z) = eps_N I(z),   photons  G(z) = (2/3) eps_pi int_z^1 I(y) dy/y,   z = E/E_X.
# The hadron spectrum is taken as x^-p (1 - x)^q. The paper gives p = 1.9 for the part that matters at
# 1e14 GeV; q sets how fast it dies toward x = 1, and p and q are varied to show what depends on them.
eN <- 0.12; ePi <- 0.73
Dh <- function(x, p, q) x^-p*(1 - x)^q
I1 <- function(x, p, q) vapply(x, function(a) integrate(function(y) Dh(y, p, q)/y, a, 1, rel.tol = 1e-10)$value, 0)
I2 <- function(x, p, q) vapply(x, function(a) integrate(function(y) I1(y, p, q)/y, a, 1, rel.tol = 1e-8)$value, 0)
above <- function(f, z) integrate(f, z, 1, rel.tol = 1e-7)$value         # number above a share z
rest <- function(x, p, q) (2/3)*ePi*I1(x, p, q)/(eN*Dh(x, p, q))
lab  <- function(z, p, q) (2/3)*ePi*I2(z, p, q)/(eN*I1(z, p, q))
labint <- function(z, p, q) (2/3)*ePi*above(function(y) I2(y, p, q), z)/(eN*above(function(y) I1(y, p, q), z))
cat("1. The quoted case: parent at rest, x = 1e-3, pure power law (the paper prints 2 to 3)\n")
for (p in c(1.5, 1.9)) cat(sprintf("   p = %.1f: photons per nucleon %.2f (q = 2), %.2f (q = 0, no cut-off); analytic (2/3)(eps_pi/eps_N)/p = %.2f\n",
                                   p, rest(1e-3, p, 2), rest(1e-3, p, 0), (2/3)*ePi/eN/p))
stopifnot(abs(rest(1e-3, 1.9, 0) - (2/3)*ePi/eN/1.9) < 0.02)
cat("\n2. A moving parent: photons per nucleon at the same energy E = z E_X\n")
cat("   z      | p=1.5 q=1  p=1.5 q=2  p=1.5 q=3 | p=1.9 q=1  p=1.9 q=2  p=1.9 q=3\n")
zs <- c(1e-3, 0.01, 0.03, 0.1, 0.2, 0.3, 0.5, 0.7)
tab <- sapply(zs, function(z) c(lab(z, 1.5, 1), lab(z, 1.5, 2), lab(z, 1.5, 3), lab(z, 1.9, 1), lab(z, 1.9, 2), lab(z, 1.9, 3)))
for (i in seq_along(zs)) cat(sprintf("   %-6g | %-9.2f  %-9.2f  %-9.2f | %-9.2f  %-9.2f  %-9.2f\n", zs[i], tab[1, i], tab[2, i], tab[3, i], tab[4, i], tab[5, i], tab[6, i]))
cat("\n3. What a search counts: photons above an energy for each nucleon above the same energy\n")
cat("   z      | p=1.5 q=2  p=1.9 q=2  p=1.9 q=3\n")
ti <- sapply(c(0.01, 0.03, 0.1, 0.2, 0.3, 0.5), function(z) c(labint(z, 1.5, 2), labint(z, 1.9, 2), labint(z, 1.9, 3)))
for (i in 1:6) cat(sprintf("   %-6g | %-9.2f  %-9.2f  %-9.2f\n", c(0.01, 0.03, 0.1, 0.2, 0.3, 0.5)[i], ti[1, i], ti[2, i], ti[3, i]))
lo <- min(ti[, 3:6]); hi <- max(ti[, 1:2])
cat(sprintf("\n   for a proton carrying a tenth to a half of the escapee's energy: %.1f to %.1f photons each; for a hundredth to a thirtieth: %.1f to %.1f\n",
            min(ti[, 3:6]), max(ti[, 3:6]), min(ti[, 1:2]), max(ti[, 1:2])))
cat(sprintf("   with eps_N at 0.10 or 0.14 and eps_pi at 0.70 or 0.76 every figure moves by a factor of %.2f to %.2f\n",
            (0.70/0.14)/(ePi/eN), (0.76/0.10)/(ePi/eN)))
# planted failure: forget the second smearing (use the at-rest spectra for a moving parent) and the
# large-z ratio must come out different, or section 2 is not testing the boost at all
d <- abs(rest(0.3, 1.9, 2)/lab(0.3, 1.9, 2) - 1)
cat(sprintf("\n   check: at z = 0.3 the at-rest and moving-parent ratios differ by %.0f per cent; at z = 1e-3 by %.1f per cent\n",
            100*d, 100*abs(rest(1e-3, 1.9, 2)/lab(1e-3, 1.9, 2) - 1)))
stopifnot(d > 0.1, abs(rest(1e-3, 1.9, 2)/lab(1e-3, 1.9, 2) - 1) < 0.05)
