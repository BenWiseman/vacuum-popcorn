# p5_observables.R: what pops would look like from Earth, if the pockets are part of the Galaxy's dark
# matter and each pop turns a fraction f of the pocket's mass-energy into chi at energy E_star.
# Inputs, all marked: E_star = 1e11 GeV (1e20 eV), the pop output; median E_out/(m_out^2/m_in) from
# p3_spectrum.R; NFW halo with r_s = 20 kpc and 0.4 GeV/cm^3 at 8.2 kpc (a common choice, not fitted
# here); Auger surface-detector photon flux limit above 40 EeV, 1.72e-4 per km^2 sr yr (scouts/A).
kpc <- 3.0857e21; yr <- 3.15576e7; GeV_g <- 1.78266e-24; tH <- 1.451e10*yr
E_star <- 1e11
cat("1. WHAT MASSES GIVE 1e20 eV (E_out = s * m_out^2/m_in)\n")
for (s in c(0.01, 0.001)) for (mo in c(1e2, 1e3, 1e6, 1e9))
  cat(sprintf("   median share s = %-6g m_out = %-8.0e GeV: m_in = %.2g GeV\n", s, mo, s*mo^2/E_star))
cat(sprintf("   and for 1e14 eV (where CHICOS searched): m_out = 100 GeV needs m_in = %.2g GeV at s = 0.01\n", 0.01*1e2^2/1e5))
cat(sprintf("   hot gas with m_out/E_typ = 1e4 and the median share 0.025 of p6: a median of 1e20 eV needs m_out = %.1e GeV (E_typ = %.0e GeV)\n",
            E_star/(0.025*1e4), E_star/(0.025*1e4)/1e4))
# halo: all-sky average of the dark-matter column, integral of rho along the line of sight to 300 kpc
rs <- 20; Rsun <- 8.2; rho_loc <- 0.4
nfw <- function(r) 1/((r/rs)*(1 + r/rs)^2)
rho_s <- rho_loc/nfw(Rsun)
col <- function(cb) integrate(function(l) rho_s*nfw(sqrt(Rsun^2 + l^2 - 2*Rsun*l*cb)), 0, 300,
                              subdivisions = 2000L, rel.tol = 1e-8)$value*kpc   # GeV/cm^2
cbs <- seq(-1, 0.999, length.out = 400)
Dbar <- mean(sapply(cbs, col))                      # uniform in cos(angle from the Galactic centre)
cat(sprintf("\n   inputs: NFW scale radius %g kpc, %g GeV/cm^3 at %g kpc; Auger photon limit above 40 EeV %.3g per km^2 sr yr (2209.05926)\n",
            rs, rho_loc, Rsun, 1.72e-4))
cat(sprintf("\n2. HALO: all-sky mean dark-matter column %.3g GeV/cm^2 (anticentre %.3g, 90 deg %.3g)\n",
            Dbar, col(-1), col(0)))
# photon cap for chi -> gamma gamma, each photon at E_star/2 = 50 EeV, above Auger's 40 EeV limit
Phi_lim <- 1.72e-4/(1e10*yr)                        # per cm^2 s sr
fG_cap <- Phi_lim*4*pi*E_star/(2*Dbar)              # per second, energy fraction per pocket per time
cat(sprintf("\n3. PHOTON CAP (chi -> 2 photons): f Gamma <= %.2g per s = %.2g per yr;\n   so at most %.1g of the pockets' mass-energy may pop per Hubble time\n",
            fG_cap, fG_cap*yr, fG_cap*tH))
# the same rate seen as an equivalent superheavy-dark-matter lifetime
cat(sprintf("   (as a decaying-dark-matter lifetime: 1/(f Gamma) = %.2g s)\n", 1/fG_cap))
# 4. one pop: how far away it can be and still put N events in a 3000 km^2 array
A <- 3000*1e10
cat("\n4. ONE POP: distance at which it leaves one event in a 3000 km^2 array (f = 1, every chi seen)\n")
for (Mg in c(1e10, 1e15, 1e20, 1e25)) {
  Nchi <- Mg/GeV_g/E_star; d1 <- sqrt(Nchi*A/(4*pi))
  cat(sprintf("   pocket %.0e g: %.2g chi, one event out to %.3g kpc\n", Mg, Nchi, d1/kpc))
}
# 5. multiplets at the photon cap, pockets = all local dark matter: rate of pops close enough for >= 1
# expected event, and the share of events that come in groups (Euclidean nearby, halo-limited far)
cat("\n5. POPS CLOSE ENOUGH TO GIVE ONE EVENT OR MORE, AT THE PHOTON CAP (all local dark matter in pockets)\n")
rho_g <- rho_loc*GeV_g                               # g/cm^3
for (Mg in c(1e10, 1e15, 1e20, 1e25)) {
  Nchi <- Mg/GeV_g/E_star; d1 <- sqrt(Nchi*A/(4*pi)); n_p <- rho_g/Mg
  R1 <- n_p*fG_cap*yr*(4*pi/3)*d1^3                 # per year, f = 1
  cat(sprintf("   pocket %.0e g: %.2g such pops per year\n", Mg, R1))
}
# 6. arrival spread of one burst at distance d: d / (2 gamma_chi^2 c), gamma_chi = E_star/m_out
cat("\n6. ARRIVAL SPREAD OF ONE BURST FROM 1 kpc\n")
for (mo in c(1e3, 1e6, 1e9)) cat(sprintf("   m_out = %.0e GeV: gamma_chi = %.1e, spread %.2g s\n", mo, E_star/mo,
                                         kpc/2.99792458e10/(2*(E_star/mo)^2)))
stopifnot(Dbar > 1e22, Dbar < 1e24)
