# p12_halo_directions.R: how many of the published cosmic rays above 100 EeV could come from pops,
# judged by arrival direction alone.
#
# If the pockets are the Galaxy's dark matter, what a pop sends arrives along the dark-matter column,
# which is several times larger toward the Galactic centre than away from it. Photons travel straight,
# and a proton of 1e20 eV bends by a few degrees in the Galactic field, far less than the width of the
# halo's pattern, so the test holds for either. Each event is drawn from
#     (1 - f) * exposure  +  f * exposure * column,
# each normalised over its own experiment's sky, and f is the share of events that follow the halo.
# Halo as in p5_observables.R: NFW, scale radius 20 kpc, Sun at 8.2 kpc, column to 300 kpc.
# Exposure: Sommers (2001) for each site, as in the void test of the parent project.
set.seed(20261002)
ev <- read.csv("../data/uhecr_above_100EeV.csv", stringsAsFactors = FALSE)
ev <- ev[!is.na(ev$ra_deg) & !is.na(ev$dec_deg), ]
cat(sprintf("%d events at or above 100 EeV with published coordinates: %s\n", nrow(ev),
            paste(names(table(ev$experiment)), table(ev$experiment), collapse = ", ")))
stopifnot(nrow(ev) == 57)
d2r <- pi/180
unitv <- function(ra, dec) cbind(cos(dec*d2r)*cos(ra*d2r), cos(dec*d2r)*sin(ra*d2r), sin(dec*d2r))
gc <- unitv(266.405, -28.936)                       # Galactic centre, J2000
rs <- 20; Rsun <- 8.2
nfw <- function(r) 1/((r/rs)*(1 + r/rs)^2)
col1 <- function(cb) integrate(function(l) nfw(sqrt(Rsun^2 + l^2 - 2*Rsun*l*cb)), 0, 300,
                               subdivisions = 4000L, rel.tol = 1e-8)$value
cbs <- c(seq(-1, 0.99, by = 0.005), 1 - 10^seq(-2, -7, length.out = 60))
colf <- splinefun(cbs, log(vapply(cbs, col1, 0)), method = "natural")
column <- function(ra, dec) exp(colf(pmin(as.vector(unitv(ra, dec) %*% t(gc)), 1 - 1e-7)))
cat(sprintf("halo column relative to the anticentre: 90 degrees from the centre %.2f, 30 degrees %.2f, 10 degrees %.2f\n",
            exp(colf(0) - colf(-1)), exp(colf(cos(30*d2r)) - colf(-1)), exp(colf(cos(10*d2r)) - colf(-1))))
expo <- function(dec, lat, thm) {
  a0 <- lat*d2r; th <- thm*d2r; d <- dec*d2r
  xi <- (cos(th) - sin(a0)*sin(d))/(cos(a0)*cos(d))
  am <- ifelse(xi > 1, 0, ifelse(xi < -1, pi, acos(pmin(1, pmax(-1, xi)))))
  pmax(cos(a0)*cos(d)*sin(am) + am*sin(a0)*sin(d), 0)
}
site <- list(TelescopeArray = c(39.3, 55), AGASA = c(35.8, 45), PierreAuger = c(-35.2, 80),
             FlysEye = c(40.2, 60), Yakutsk = c(61.7, 60))
n <- 40000; k <- 0:(n - 1); zz <- 1 - 2*(k + 0.5)/n; phi <- k*pi*(3 - sqrt(5))
g_ra <- (phi/d2r) %% 360; g_dec <- asin(zz)/d2r; g_col <- column(g_ra, g_dec)
# per experiment: mean of the column over its exposure, so that column/mean integrates to one
norm <- lapply(site, function(s) { w <- expo(g_dec, s[1], s[2]); list(w = w, mean_col = sum(w*g_col)/sum(w)) })
ratio <- function(e) column(e$ra_deg, e$dec_deg)/vapply(e$experiment, function(x) norm[[x]]$mean_col, 0)
loglik <- function(f, r) sum(log(1 - f + f*r))
fit <- function(r) {
  o <- optimize(function(f) -loglik(f, r), c(0, 1)); fh <- if (loglik(0, r) >= -o$objective) 0 else o$minimum
  up <- tryCatch(uniroot(function(f) 2*(loglik(fh, r) - loglik(f, r)) - 2.71, c(fh, 1))$root, error = function(e) 1)
  c(f = fh, up = up, ts = 2*(loglik(fh, r) - loglik(0, r)))
}
draw <- function(m, exper, f) {                     # m directions for one experiment, a share f from the halo
  w <- norm[[exper]]$w; nh <- rbinom(1, m, f)
  i <- c(sample(n, nh, replace = TRUE, prob = w*g_col), sample(n, m - nh, replace = TRUE, prob = w))
  data.frame(experiment = exper, ra_deg = g_ra[i], dec_deg = g_dec[i])
}
fake <- function(f) do.call(rbind, lapply(names(table(ev$experiment)), function(x) draw(sum(ev$experiment == x), x, f)))

cat("\nVALIDATION (fake skies with the real number of events per experiment)\n")
null <- t(replicate(400, fit(ratio(fake(0)))))
cat(sprintf("   no halo events: fitted share is 0 in %.0f per cent of skies; the 95 per cent upper limit lies below the truth (0) in none, and has median %.2f\n",
            100*mean(null[, "f"] == 0), median(null[, "up"])))
for (ft in c(0.2, 0.5)) { p <- t(replicate(400, fit(ratio(fake(ft)))))
  cat(sprintf("   planted share %.1f: median fitted %.2f; upper limit below the truth in %.1f per cent of skies (5 expected)\n",
              ft, median(p[, "f"]), 100*mean(p[, "up"] < ft)))
  stopifnot(abs(median(p[, "f"]) - ft) < 0.08, mean(p[, "up"] < ft) < 0.10) }
p1 <- t(replicate(400, fit(ratio(fake(1)))))
cat(sprintf("   every event from the halo: median fitted share %.2f, and the test statistic beats the largest of the 400 no-halo skies in %.0f per cent\n",
            median(p1[, "f"]), 100*mean(p1[, "ts"] > max(null[, "ts"]))))
stopifnot(median(p1[, "f"]) > 0.9)

cat("\nDATA\n")
r <- ratio(ev); a <- fit(r)
cat(sprintf("   all %d events: fitted halo share %.2f, 95 per cent upper limit %.2f, that is %.0f events; test statistic %.2f (no-halo skies exceed it %.0f per cent of the time)\n",
            nrow(ev), a["f"], a["up"], a["up"]*nrow(ev), a["ts"], 100*mean(null[, "ts"] >= a["ts"])))
for (x in c("PierreAuger", "TelescopeArray", "AGASA")) { j <- ev$experiment == x; b <- fit(r[j])
  cat(sprintf("   %-14s %2d events: fitted share %.2f, upper limit %.2f (%.0f events); mean column over the sky's mean %.2f\n",
              x, sum(j), b["f"], b["up"], b["up"]*sum(j), mean(r[j]))) }
big <- ev[ev$energy_EeV >= 200, ]
cat("   the named events with coordinates, column relative to their experiment's sky average:\n")
for (i in seq_len(nrow(big))) cat(sprintf("      %s %s %g EeV: %.2f, %.0f degrees from the Galactic centre\n", big$experiment[i], big$date[i],
    big$energy_EeV[i], ratio(big[i, ]), acos(unitv(big$ra_deg[i], big$dec_deg[i]) %*% t(gc))/d2r))
cat(sprintf("   Auger events within 30 degrees of the Galactic centre: %d of %d; expected from exposure alone %.1f, if all came from the halo %.1f\n",
    sum(acos(unitv(ev$ra_deg, ev$dec_deg) %*% t(gc))[ev$experiment == "PierreAuger"]/d2r < 30), sum(ev$experiment == "PierreAuger"),
    35*sum(norm$PierreAuger$w*(acos(unitv(g_ra, g_dec) %*% t(gc))/d2r < 30))/sum(norm$PierreAuger$w),
    35*sum((norm$PierreAuger$w*g_col)*(acos(unitv(g_ra, g_dec) %*% t(gc))/d2r < 30))/sum(norm$PierreAuger$w*g_col)))

cat("\nBY EXPERIMENT, against fake skies drawn from that experiment's exposure alone\n")
for (x in c("PierreAuger", "TelescopeArray", "AGASA")) { j <- ev$experiment == x; m <- sum(j)
  sim <- replicate(4000, mean(ratio(draw(m, x, 0)))); simh <- replicate(2000, mean(ratio(draw(m, x, 1))))
  cat(sprintf("   %-14s mean column over sky average %.2f; no-halo skies give %.2f +- %.2f and reach the observed value %.1f per cent of the time; all-halo skies give %.2f +- %.2f and fall to it %.1f per cent of the time\n",
              x, mean(r[j]), mean(sim), sd(sim), 100*mean(sim >= mean(r[j])), mean(simh), sd(simh), 100*mean(simh <= mean(r[j])))) }
ta <- ev[ev$experiment == "TelescopeArray", ]
cat("   Telescope Array events (date, EeV, degrees from the Galactic centre, column over sky average):\n")
for (i in order(-ratio(ta))) cat(sprintf("      %s %5.0f %4.0f %.2f\n", ta$date[i], ta$energy_EeV[i],
    acos(unitv(ta$ra_deg[i], ta$dec_deg[i]) %*% t(gc))/d2r, ratio(ta[i, ])))

cat("\nCONTROL: the events between 57 and 100 EeV, same test\n")
ct <- read.csv("../data/uhecr_57_to_100EeV_control.csv", stringsAsFactors = FALSE)
ct <- ct[!is.na(ct$ra_deg) & !is.na(ct$dec_deg) & ct$experiment %in% names(site), ]
rc <- ratio(ct)
for (x in intersect(c("PierreAuger", "TelescopeArray", "AGASA"), unique(ct$experiment))) { j <- ct$experiment == x; m <- sum(j)
  sim <- replicate(2000, mean(ratio(draw(m, x, 0))))
  cat(sprintf("   %-14s %3d events: mean column over sky average %.2f; no-halo skies give %.2f +- %.2f and reach it %.1f per cent of the time; fitted share %.2f, upper limit %.2f\n",
              x, m, mean(rc[j]), mean(sim), sd(sim), 100*mean(sim >= mean(rc[j])), fit(rc[j])["f"], fit(rc[j])["up"])) }
cat("\nENERGY: Telescope Array events above 57 EeV, column against energy\n")
tall <- rbind(ev[ev$experiment == "TelescopeArray", c("energy_EeV", "ra_deg", "dec_deg", "experiment", "date")],
              ct[ct$experiment == "TelescopeArray", c("energy_EeV", "ra_deg", "dec_deg", "experiment", "date")])
tall <- tall[tall$date < "2014-01-01", ]            # the complete 2008-2013 sample, without the later Amaterasu
rt <- ratio(tall); sp <- cor(tall$energy_EeV, rt, method = "spearman")
perm <- replicate(20000, cor(sample(tall$energy_EeV), rt, method = "spearman"))
cat(sprintf("   %d events of 2008 to 2013: Spearman correlation of energy with halo column %.2f; shuffled energies reach it %.2f per cent of the time\n",
            nrow(tall), sp, 100*mean(perm >= sp)))
for (cut in c(57, 70, 80, 100, 120)) { j <- tall$energy_EeV >= cut
  cat(sprintf("   at or above %3d EeV: %2d events, mean column over sky average %.2f, within 50 degrees of the Galactic centre %d\n", cut, sum(j), mean(rt[j]),
      sum(acos(unitv(tall$ra_deg[j], tall$dec_deg[j]) %*% t(gc))/d2r < 50))) }
frac50 <- sum(norm$TelescopeArray$w*(acos(unitv(g_ra, g_dec) %*% t(gc))/d2r < 50))/sum(norm$TelescopeArray$w)
cat(sprintf("   share of Telescope Array's exposure within 50 degrees of the Galactic centre: %.3f\n", frac50))

cat("\nWHAT WOULD SETTLE THE LEAN: Telescope Array's events above 100 EeV with no published direction\n")
ang_g <- acos(pmin(1, unitv(g_ra, g_dec) %*% t(gc)))/d2r; wt <- norm$TelescopeArray$w
p_iso <- sum(wt*(ang_g < 50))/sum(wt); p_halo <- sum(wt*g_col*(ang_g < 50))/sum(wt*g_col)
t10 <- ev[ev$experiment == "TelescopeArray" & ev$date < "2014-01-01", ]
a10 <- acos(unitv(t10$ra_deg, t10$dec_deg) %*% t(gc))/d2r; k10 <- sum(a10 < 50)
sim10 <- replicate(20000, mean(ratio(draw(10, "TelescopeArray", 0))))
cat(sprintf("   the complete 2008 to 2013 sample: %d events, %d within 50 degrees of the Galactic centre; mean column over sky average %.2f, which no-halo skies reach %.2f per cent of the time\n",
            nrow(t10), k10, mean(ratio(t10)), 100*mean(sim10 >= mean(ratio(t10)))))
cat(sprintf("   chance that a Telescope Array event lies within 50 degrees of the centre: %.3f if isotropic, %.3f if every event followed the halo; chance of %d or more in %d if isotropic %.4f\n",
            p_iso, p_halo, k10, nrow(t10), 1 - pbinom(k10 - 1, nrow(t10), p_iso)))
more <- 28 - 11                                      # 2311.14231 counts 28 above 100 EeV; 11 have published directions
cat(sprintf("   of the %d further events: %.1f expected within 50 degrees if isotropic, %.1f if every event followed the halo, %.1f at the rate of the 2008 to 2013 sample\n",
            more, more*p_iso, more*p_halo, more*k10/nrow(t10)))
cat(sprintf("   chance of 4 or more among %d if isotropic: %.3f; of 1 or fewer at the published rate: %.3f\n",
            more, 1 - pbinom(3, more, p_iso), pbinom(1, more, k10/nrow(t10))))
