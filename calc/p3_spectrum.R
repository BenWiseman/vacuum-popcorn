# p3_spectrum.R: what a collapsing pocket throws out, for particles initially at rest and spread
# uniformly through it. Energies in units of m_in; mu = m_out/m_in; results quoted as E/(m_out^2/m_in).
# xmin = (wall thickness)/R0 is how far the pocket shrinks before its wall unravels; a particle still
# inside then is counted as "trapped at the end" (its fate is set by the wall's final collapse, which
# this calculation does not follow).
source("pt_lib.R")
args <- commandArgs(trailingOnly = TRUE)
N <- if (length(args)) as.integer(args[1]) else 1500
set.seed(17)
q <- function(v, p) if (length(v)) quantile(v, p, names = FALSE) else NA
run <- function(mu, k, xmin, N) {
  x0 <- runif(N)^(1/3); x0 <- pmin(x0, 1 - 1e-9)
  res <- lapply(x0, function(x) track(x, mu, k, xmin = xmin, max_events = 2000))
  st <- vapply(res, `[[`, "", "status"); E <- vapply(res, `[[`, 0, "E")
  first <- vapply(res, function(r) as.character(r$first), "")
  nc <- vapply(res, `[[`, 0, "ncatch")
  esc <- st == "escaped"; s <- E[esc]/mu^2
  list(mu = mu, k = k, xmin = xmin, n = N, f_esc = mean(esc), f_trap = mean(st == "trapped_at_end"),
       f_other = mean(!(st %in% c("escaped", "trapped_at_end"))),
       f_swept_first = mean(first == "transmit"),
       q10 = q(s, 0.1), q50 = q(s, 0.5), q90 = q(s, 0.9), smax = if (any(esc)) max(s) else NA,
       ceiling = (2*mu^2 - 1)/mu^2, emean = if (any(esc)) mean(s) else NA,
       f_above = mean(s > (2*mu^2 - 1)/mu^2), f_catch = mean(nc[esc] > 0),
       s_all = s, x0 = x0, status = st, ncatch = nc)
}
out <- list()
cat(sprintf("%-6s %-5s %-7s | %-6s %-6s %-6s %-6s %-6s | %-9s %-9s %-9s %-8s %-8s %-7s\n", "mu", "k", "xmin",
            "esc", "trap", "other", "swept1", "caught", "E10%", "E50%", "E90%", "Emean", "Emax", "above2"))
for (mu in c(1e2, 1e4, 1e6, 1e9)) for (k in c(0.1, 1, 10)) for (xmin in c(1e-6, 1e-12, 1e-20)) {
  r <- run(mu, k, xmin, N); out[[length(out) + 1]] <- r
  cat(sprintf("%-6.0e %-5g %-7.0e | %-6.3f %-6.3f %-6.3f %-6.3f %-6.3f | %-9.3g %-9.3g %-9.3g %-8.3g %-8.3g %-7.4f\n",
              mu, k, xmin, r$f_esc, r$f_trap, r$f_other, r$f_swept_first, r$f_catch, r$q10, r$q50, r$q90,
              r$emean, r$smax, r$f_above))
  flush.console()
}
saveRDS(out, "p3_spectrum.rds")
cat("\nEnergies are E_out / (m_out^2 / m_in) for escaped particles. above2 = share above the head-on\nceiling (2 m_out^2 - m_in^2)/m_in; caught = share of escapers caught from behind at least once.\n")
