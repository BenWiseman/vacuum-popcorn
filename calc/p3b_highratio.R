# p3b_highratio.R: the gas-at-rest spectrum at the high mass ratios the letter's examples need
# (a blind judge noticed that "1 TeV outside, 100 eV inside" is a ratio of 1e10, beyond p3's 1e9).
# Also: how small x_min = (wall thickness)/R0 is for the pockets of section 5, since the median
# depends on it at high ratios.
source("pt_lib.R")
set.seed(53); N <- 1000
q <- function(v, p) if (length(v)) quantile(v, p, names = FALSE) else NA
cat(sprintf("%-6s %-5s %-7s | %-6s %-6s | %-9s %-9s %-9s\n", "mu", "k", "xmin", "esc", "trap", "E10%", "E50%", "E90%"))
for (mu in c(1e6, 1e8, 1e9, 1e10, 1e11, 1e12)) for (k in c(0.1, 1, 10)) for (xmin in c(1e-12, 1e-15)) {
  x0 <- pmin(runif(N)^(1/3), 1 - 1e-9)
  res <- lapply(x0, function(x) track(x, mu, k, xmin = xmin, max_events = 2000))
  st <- vapply(res, `[[`, "", "status"); E <- vapply(res, `[[`, 0, "E")
  s <- E[st == "escaped"]/mu^2
  cat(sprintf("%-6.0e %-5g %-7.0e | %-6.3f %-6.3f | %-9.3g %-9.3g %-9.3g\n", mu, k, xmin, mean(st == "escaped"),
              mean(st == "trapped_at_end"), q(s, 0.1), q(s, 0.5), q(s, 0.9)))
  flush.console()
}
# x_min for the pockets of section 5: wall thickness ~ 1/L, radius from the mass at energy density L^4
GeV_g <- 1.78266e-24
cat("\nx_min = 1/(L R0) for a pocket of mass M and vacuum scale L (R0 from M = (4 pi/3) L^4 R0^3)\n")
for (Mg in c(1e15, 1e20)) for (L in c(1e-3, 1, 1e3)) {
  M <- Mg/GeV_g; R0 <- (3*M/(4*pi*L^4))^(1/3)
  cat(sprintf("   M = %.0e g, L = %g GeV: x_min = %.1e\n", Mg, L, 1/(L*R0)))
}
