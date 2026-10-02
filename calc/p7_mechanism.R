# p7_mechanism.R: is the bulk of the spectrum the two-step path the findings describe? A particle at
# rest is reflected once by a still-slow wall (energy E1), races through the centre, and crosses the
# far side, where the wall is ultra-fast, with E1 + m_out^2/(4 E1) to m_out^2/(2 E1) (the head-on
# fast-wall result, which assumes the far-side wall is ultra-fast). Units m_in = 1.
source("pt_lib.R")
set.seed(31)
for (cfg in list(c(1e6, 1), c(1e6, 0.1), c(1e4, 10))) {
  mu <- cfg[1]; k <- cfg[2]; N <- 1500
  x0 <- pmin(runif(N)^(1/3), 1 - 1e-9)
  res <- lapply(x0, function(x) track(x, mu, k, xmin = 1e-12, max_events = 2000))
  st <- vapply(res, `[[`, "", "status"); n <- vapply(res, `[[`, 0, "n"); E <- vapply(res, `[[`, 0, "E")
  E1 <- vapply(res, `[[`, 0, "E1"); nc <- vapply(res, `[[`, 0, "ncatch")
  two <- st == "escaped" & n == 2 & nc == 0
  lo <- E1 + mu^2/(4*E1); hi <- mu^2/(2*E1)
  inband <- two & E >= lo*(1 - 1e-9) & E <= hi*(1 + 1e-9)
  cat(sprintf("mu %.0e k %-4g: escaped %.3f; two-step path (reflect, then cross the far side) %.3f of escapers;\n   of those, inside [E1 + m_out^2/(4E1), m_out^2/(2E1)]: %.3f; their share of the escaped energy %.3f\n",
              mu, k, mean(st == "escaped"), sum(two)/sum(st == "escaped"), sum(inband)/max(sum(two), 1),
              sum(E[two])/sum(E[st == "escaped"])))
  stopifnot(sum(inband) >= 0.85*sum(two))   # the band assumes the far-side wall is ultra-fast; a slower one gives a little less
}
