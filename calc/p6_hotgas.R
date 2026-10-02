# p6_hotgas.R: what a collapsing pocket throws when its trapped gas is hot or degenerate.
# Isotropic momenta, positions uniform in the pocket. Units m_in = 1. Two gases:
#   thermal, relativistic: |p| drawn from p^2 exp(-p/T) (the massless limit, T >> m_in);
#   degenerate: |p| uniform in the Fermi ball, |p| = p_F u^(1/3).
# Results quoted as E_out / (m_out^2 / E_typ), with E_typ = 3T (thermal mean) or 3 p_F / 4 (Fermi mean).
source("pt3_lib.R")
args <- commandArgs(trailingOnly = TRUE); N <- if (length(args)) as.integer(args[1]) else 400
set.seed(41)
q <- function(v, pr) if (length(v)) quantile(v, pr, names = FALSE) else NA
cat(sprintf("%-8s %-7s %-7s %-4s | %-6s %-6s %-6s | %-9s %-9s %-9s %-9s\n", "gas", "E_typ", "mout/E", "k",
            "esc", "trap", "caught", "E10%", "E50%", "E90%", "Emax"))
tabs <- lapply(c(0.1, 1, 10), make_lagtab)
for (gas in c("thermal", "fermi")) for (Et in c(10, 1000)) for (ratio in c(1e2, 1e4)) for (ti in 1:3) {
  tab <- tabs[[ti]]; mu <- ratio*Et
  res <- lapply(1:N, function(i) {
    x0 <- runif(1)^(1/3)*(1 - 1e-9); d <- rnorm(3); r0 <- x0*d/sqrt(sum(d^2))
    pm <- if (gas == "thermal") rgamma(1, shape = 3, scale = Et/3) else (4*Et/3)*runif(1)^(1/3)
    w <- rnorm(3); p0 <- pm*w/sqrt(sum(w^2))
    track3(r0, p0, mu, tab, xmin = 1e-12)
  })
  st <- vapply(res, `[[`, "", "status"); E <- vapply(res, `[[`, 0, "E"); nc <- vapply(res, `[[`, 0, "ncatch")
  esc <- st == "escaped"; s <- E[esc]/(mu^2/Et)
  cat(sprintf("%-8s %-7g %-7.0e %-4g | %-6.3f %-6.3f %-6.3f | %-9.3g %-9.3g %-9.3g %-9.3g\n", gas, Et, ratio, tab$k,
              mean(esc), mean(st == "trapped_at_end"), mean(nc[esc] > 0), q(s, 0.1), q(s, 0.5), q(s, 0.9),
              if (any(esc)) max(s) else NA))
}
cat("\nEnergies are E_out / (m_out^2 / E_typ) for escaped particles.\n")
