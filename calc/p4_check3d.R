# p4_check3d.R: the 3D tracker against the radial one, for particles at rest (where both apply).
source("pt3_lib.R")
set.seed(23); worst <- 0; agree <- 0; N <- 0
for (k in c(0.1, 1, 10)) {
  tab <- make_lagtab(k)
  for (mu in c(1e2, 1e4, 1e6)) for (x0 in c(0.2, 0.5, 0.8, 0.95)) {
    dirn <- rnorm(3); dirn <- dirn/sqrt(sum(dirn^2))
    a <- track(x0, mu, k, xmin = 1e-12, max_events = 2000)
    b <- track3(x0*dirn, c(0, 0, 0), mu, tab, xmin = 1e-12)
    N <- N + 1; same <- a$status == b$status; if (same) agree <- agree + 1
    rel <- if (same && a$status == "escaped") abs(a$E/b$E - 1) else NA
    if (!is.na(rel)) worst <- max(worst, rel)
    cat(sprintf("k %4g mu %.0e x0 %.2f: radial %-15s %-11.6g | 3D %-15s %-11.6g %s\n", k, mu, x0,
                a$status, a$E, b$status, b$E, if (is.na(rel)) "" else sprintf("rel %.1e", rel)))
  }
}
cat(sprintf("\n%d of %d agree on the outcome; largest relative energy difference %.1e\n", agree, N, worst))
stopifnot(agree == N, worst < 1e-6)
