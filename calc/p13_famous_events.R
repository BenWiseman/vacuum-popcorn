# p13_famous_events.R: how many of the events above 100 EeV could pops account for, given each
# experiment's own photon search. Inputs are literature values from scouts/O_famous_events_primary_type.md
# and scouts/A_amaterasu_and_photon_limits.md (read by a scout, not re-opened here), each named below.
#
# A pop's escapees decay; if they decay to quarks the jets hold photons and nucleons, and every pop
# proton (or neutron: from the halo a neutron of 1e20 eV arrives before it decays) comes with photons.
# How many is worked out in p14_photons_per_proton.R from the fragmentation ratios of hep-ph/0307279:
# about 2.2 for each nucleon when the nucleon carries a hundredth of the escapee's energy, falling to
# about 0.5 when it carries half. (The "2 to 3" usually quoted is for a particle decaying at rest and
# a product carrying a thousandth of its energy; an earlier version of this script used it.) A photon
# search that found none then caps the protons.
gN <- c(0.5, 2.2)
# Telescope Array: photon flux limit above 1e20 eV 3.0e-4 per km2 sr yr with 0 candidates, effective
# photon exposure 10.44e3 km2 sr yr (2512.01638); total exposure 1.6e4 with 28 events above 100 EeV
# (2311.14231).
ta_lim <- 3.0e-4; ta_phot <- 10.44e3; ta_tot <- 1.6e4; ta_n <- 28
# Auger: photon flux limit above 40 EeV 1.72e-4 with 0 candidates (2209.05926); total exposure 122,000
# (2206.13492); 35 events at or above 100 EeV in the Phase I catalog (data/uhecr_above_100EeV.csv).
au_lim <- 1.72e-4; au_tot <- 122000; au_n <- 35
# AGASA and Yakutsk: photon fraction below 0.36 above 1e20 eV, 95 per cent (astro-ph/0601449); 11 AGASA events.
ag_frac <- 0.36; ag_n <- 11
cat(sprintf("1. Most pop protons each photon null allows (95 per cent), for %.1f and for %.1f photons per proton\n", gN[2], gN[1]))
cat(sprintf("   Telescope Array: %.1f to %.1f of its %d events above 100 EeV\n", ta_lim/gN[2]*ta_tot, ta_lim/gN[1]*ta_tot, ta_n))
cat(sprintf("   Auger:           %.1f to %.1f of its %d\n", au_lim/gN[2]*au_tot, au_lim/gN[1]*au_tot, au_n))
cat(sprintf("   AGASA:           %.1f to %.1f of its %d\n", ag_frac/gN[2]*ag_n, ag_frac/gN[1]*ag_n, ag_n))
cat(sprintf("   together:        %.0f to %.0f of %d\n", ta_lim/gN[2]*ta_tot + au_lim/gN[2]*au_tot + ag_frac/gN[2]*ag_n,
            ta_lim/gN[1]*ta_tot + au_lim/gN[1]*au_tot + ag_frac/gN[1]*ag_n, ta_n + au_n + ag_n))
cat("\n2. If k of Telescope Array's events above 100 EeV were pop protons: photons its search should have held, and the chance of the none it found\n")
for (k in c(1, 2, 4)) { mu <- k/ta_tot*gN*ta_phot
  cat(sprintf("   k = %d: %.1f to %.1f photons expected; chance of none %.3f to %.4f\n", k, mu[1], mu[2], exp(-mu[1]), exp(-mu[2]))) }
cat(sprintf("   one such proton in Telescope Array's exposure means a photon flux above 1e20 eV of %.2g to %.2g per km2 sr yr; Auger's limit above 40 EeV is %.2g\n",
            gN[1]/ta_tot, gN[2]/ta_tot, au_lim))
cat("\n3. The named events as photons (probability of a photon primary, or the printed exclusion)\n")
cat("   Fly's Eye 1991, 320 EeV:  0.13, not excluded; the profile lies 1.5 sigma from a photon's (astro-ph/0410739)\n")
cat("   Telescope Array 2021, 244 EeV: excluded at 99.986 per cent (2311.14231)\n")
cat("   AGASA 1993, 213 EeV:      0.001; none of the 10 highest-energy AGASA and Yakutsk events with muon data was a photon, at 95 per cent (astro-ph/0601449)\n")
cat("   Yakutsk 1989:             0.000 (astro-ph/0601449)\n")
stopifnot(abs(ta_lim/gN[2]*ta_tot - 3.0e-4/2.2*1.6e4) < 1e-9, exp(-1/ta_tot*gN[1]*ta_phot) > exp(-1/ta_tot*gN[2]*ta_phot))
cat("\n4. What would rule out the 244 EeV event as a pop proton: Telescope Array's photon exposure above 1e20 eV at which finding no photon has a chance below 5 per cent\n")
need <- -log(0.05)*ta_tot/rev(gN)
cat(sprintf("   %.0f to %.0f km2 sr yr, which is %.1f to %.1f times today's %.0f\n", need[1], need[2], need[1]/ta_phot, need[2]/ta_phot, ta_phot))
