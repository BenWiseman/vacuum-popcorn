# p3c_examples.R: the letter's mass examples, simulated at their own mass ratios rather than
# extrapolated. Gas at rest, pocket ending at x_min = 1e-13 (the middle of the 2e-12 to 4e-16 range
# that p3b computes for 1e15 to 1e20 g pockets), three wall histories.
source("pt_lib.R")
set.seed(61); N <- 1500
ex <- list(c(1e7, 10), c(1e3, 1e-8), c(1e3, 1e-7))     # (m_out, m_in) in GeV
for (e in ex) {
  mo <- e[1]; mi <- e[2]; mu <- mo/mi
  cat(sprintf("m_out = %.0e GeV, m_in = %.0e GeV (ratio %.0e), m_out^2/m_in = %.1e eV\n", mo, mi, mu, mo^2/mi*1e9))
  for (k in c(0.1, 1, 10)) {
    x0 <- pmin(runif(N)^(1/3), 1 - 1e-9)
    res <- lapply(x0, function(x) track(x, mu, k, xmin = 1e-13, max_events = 2000))
    st <- vapply(res, `[[`, "", "status"); E <- vapply(res, `[[`, 0, "E")*mi*1e9   # eV
    Ee <- E[st == "escaped"]
    cat(sprintf("   k = %-4g escaped %.2f | median %.1e eV, top tenth above %.1e eV, top hundredth above %.1e eV | share above 1e20 eV %.3f\n",
                k, mean(st == "escaped"), median(Ee), quantile(Ee, 0.9), quantile(Ee, 0.99), mean(Ee > 1e20)))
    flush.console()
  }
}
