# p11_loaded_summary.R: what the loaded-collapse scan (p10_loaded.R, read from p10_loaded.rds) says,
# in the terms the findings quote. Runs in a second; rerun p10 first if pl_lib.R has changed.
r <- readRDS("p10_loaded.rds")
d <- do.call(rbind, lapply(r, function(a) as.data.frame(a)))
d$run <- c(rep("A", 36), rep("B", 6), rep(c("C standard", "C fine"), 4))
stopifnot(nrow(d) == 50, all(d$mu[1:36] %in% c(1e4, 1e6)), all(d$xmin[37:42] == 1e-12), all(d$per_decade[c(44, 46, 48, 50)] == 240))
A <- d[d$run == "A", ]
cat("1. Below the load at which the wall stalls, the light-load results stand.\n")
cat("   mu     k    | heaviest load with no stall | escaped share there / at no load | median energy there / at no load | escapees' share of the pocket's energy there | particles' share\n")
for (mu in c(1e4, 1e6)) for (k in c(0.1, 1, 10)) {
  s <- A[A$mu == mu & A$k == k, ]; s0 <- s[s$ell == 0, ]; ok <- s[s$x_stall == 0 & s$ell > 0, ]; top <- ok[which.max(ok$ell), ]
  first <- s[s$x_stall > 0, ]; first <- if (nrow(first)) min(first$ell) else NA
  cat(sprintf("   %-6.0e %-4g | %3g (stalls at %s) | %.3f / %.3f | %.4f / %.4f | %.2f | %.2f\n", mu, k, top$ell,
              if (is.na(first)) "none up to 80" else as.character(first), top$f_esc, s0$f_esc, top$q50, s0$q50, top$eff, top$drained))
}
u <- A[A$x_stall == 0 & A$ell > 0, ]; base <- A[match(paste(u$mu, u$k), paste(A$mu, A$k)), ]      # each config's no-load row comes first
stopifnot(all(base$ell == 0))
cat(sprintf("   over the %d loaded collapses that never stall: escaped share changes by at most %.3f, median energy by at most %.1f per cent, 90th percentile by at most %.1f per cent\n",
            nrow(u), max(abs(u$f_esc - base$f_esc)), 100*max(abs(u$q50/base$q50 - 1)), 100*max(abs(u$q90/base$q90 - 1))))
big <- u$mu == 1e6
cat(sprintf("   of those, at mass ratio 1e4: at most %.3f, %.1f per cent, %.1f per cent; at 1e6: %.3f, %.1f per cent, %.1f per cent\n",
            max(abs(u$f_esc - base$f_esc)[!big]), 100*max(abs(u$q50/base$q50 - 1)[!big]), 100*max(abs(u$q90/base$q90 - 1)[!big]),
            max(abs(u$f_esc - base$f_esc)[big]), 100*max(abs(u$q50/base$q50 - 1)[big]), 100*max(abs(u$q90/base$q90 - 1)[big])))
cat(sprintf("   the particles' share of the pocket's energy in the collapses that never stall: up to %.2f; in those that stall: %.3f to %.3f\n",
            max(u$drained), min(A$drained[A$x_stall > 0]), max(A$drained[A$x_stall > 0])))

cat("\n2. Past it, the wall stops early and fewer get out.\n")
cat("   mu     k    load | pocket size when the wall stops | escaped share (no load) | median energy (no load) | escapees' share of the pocket's energy\n")
st <- A[A$x_stall > 0, ]
for (i in seq_len(nrow(st))) { s0 <- A[A$mu == st$mu[i] & A$k == st$k[i] & A$ell == 0, ]
  cat(sprintf("   %-6.0e %-4g %-4g | %.1e | %.3f (%.3f) | %.4f (%.4f) | %.2f\n", st$mu[i], st$k[i], st$ell[i], st$x_stall[i], st$f_esc[i], s0$f_esc, st$q50[i], s0$q50, st$eff[i])) }
cat(sprintf("   escapees' share of the pocket's energy, over every loaded collapse: largest %.2f (mu %.0e, k %g, load %g)\n",
            max(A$eff), A$mu[which.max(A$eff)], A$k[which.max(A$eff)], A$ell[which.max(A$eff)]))

cat("\n3. How small a trapped gas loads a pocket: its rest mass as a share of the pocket's energy is load/mu^2.\n")
for (mu in c(1e4, 1e6)) for (ell in c(10, 40)) cat(sprintf("   mu %.0e, load %g: %.0e\n", mu, ell, ell/mu^2))

cat("\n4. The rare escapees above m_out^2/m_in (loaded runs, where the shells are fine enough to count them).\n")
t <- A[A$ell == 3, ]
for (i in seq_len(nrow(t))) cat(sprintf("   mu %.0e k %-4g: %.2f per cent of escapees, carrying %.0f per cent of the escapees' energy; 99th percentile %.3f\n",
                                       t$mu[i], t$k[i], 100*t$f_tail[i], 100*t$e_tail[i], t$q99[i]))

cat("\n5. A deeper pocket (ending at 1e-12 of its size instead of 1e-6), mu 1e4, k 1: same loads, side by side.\n")
B <- d[d$run == "B", ]; A1 <- A[A$mu == 1e4 & A$k == 1, ]
for (i in seq_len(nrow(B))) cat(sprintf("   load %-3g: escaped %.3f (1e-6: %.3f) | median %.4f (%.4f) | escapees' share %.2f (%.2f) | wall stops at %.1e (%.1e)\n",
  B$ell[i], B$f_esc[i], A1$f_esc[i], B$q50[i], A1$q50[i], B$eff[i], A1$eff[i], B$x_stall[i], A1$x_stall[i]))

cat("\n6. Sampling: the standard run against the finer one (mu 1e4, k 1), load by load.\n")
S <- d[d$run == "C standard", ]; F <- d[d$run == "C fine", ]
for (i in 1:4) cat(sprintf("   load %-3g: escaped %.3f / %.3f | median %.4f / %.4f | escapees' share %.3f / %.3f | particles' share %.3f / %.3f | wall stops at %.1e / %.1e\n",
  S$ell[i], S$f_esc[i], F$f_esc[i], S$q50[i], F$q50[i], S$eff[i], F$eff[i], S$drained[i], F$drained[i], S$x_stall[i], F$x_stall[i]))
ns <- S$x_stall == 0 & F$x_stall == 0
cat(sprintf("   where neither stalls (load %s): the finer run gives the particles %.0f per cent more of the pocket's energy\n",
            paste(S$ell[ns], collapse = ", "), 100*max(F$drained[ns]/S$drained[ns] - 1)))
flip <- (S$x_stall == 0) != (F$x_stall == 0)
cat(sprintf("   the two disagree on whether the wall stalls at load %s; elsewhere the escaped share differs by at most %.3f\n",
            paste(S$ell[flip], collapse = ", "), max(abs(S$f_esc - F$f_esc)[!flip])))
cat(sprintf("   load at which the particles would take all the pocket's energy, from the heaviest load with no stall in each: standard %.0f, finer %.0f\n",
            max(S$ell[S$x_stall == 0])/S$drained[which.max(S$ell*(S$x_stall == 0))], max(F$ell[F$x_stall == 0])/F$drained[which.max(F$ell*(F$x_stall == 0))]))
