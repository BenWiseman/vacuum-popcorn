# p4b_brute3d.R: the 3D tracker against brute-force time stepping for particles with real momenta,
# where double precision is safe (small mu, pocket ending at 2% of its radius, so gamma stays below ~1e4).
source("pt3_lib.R")
brute3 <- function(r0, p0, mu, k, xmin = 0.02, dt0 = 5e-5) {
  X <- 1 - 1e-12; r <- r0; p <- p0; E <- sqrt(1 + sum(p^2)); m2 <- 1   # a wall exactly at rest never starts moving here
  wv <- function(X) beta_x(X, k)
  repeat {
    dt <- min(dt0, 0.001*X); v <- p/E
    Xn <- X - dt*wv(X - 0.5*dt*wv(X)); rn <- r + v*dt
    if (Xn < xmin) return(list(status = "trapped_at_end", E = E))
    if (sqrt(sum(rn^2)) >= Xn) {
      lo <- 0; hi <- dt
      for (it in 1:70) { mid <- (lo + hi)/2
        Xm <- X - mid*wv(X - 0.5*mid*wv(X)); rm <- r + v*mid
        if (sqrt(sum(rm^2)) >= Xm) hi <- mid else lo <- mid }
      Xh <- X - hi*wv(X - 0.5*hi*wv(X)); rh <- r + v*hi; nrm <- rh/sqrt(sum(rh^2))
      pn <- sum(p*nrm); ptv <- p - pn*nrm; pt2 <- sum(ptv^2)
      e <- encounter(E, pn, pt2, 1, mu, gam(Xh, k))
      if (e$type == "transmit") return(list(status = "escaped", E = e$E))
      if (e$type == "reflect") { E <- e$E; p <- e$pn*nrm + ptv }
      X <- Xh; r <- nrm*Xh*(1 - 1e-12); next
    }
    X <- Xn; r <- rn
  }
}
set.seed(29); worst <- 0; agree <- 0; N <- 0
for (k in c(0.3, 3)) {
  tab <- make_lagtab(k)
  for (mu in c(2, 4)) for (i in 1:6) {
    x0 <- runif(1)^(1/3)*0.97; d <- rnorm(3); r0 <- x0*d/sqrt(sum(d^2))
    q <- rnorm(3); p0 <- q/sqrt(sum(q^2))*10^runif(1, -1, 0.5)
    a <- track3(r0, p0, mu, tab, xmin = 0.02); b <- brute3(r0, p0, mu, k, xmin = 0.02)
    N <- N + 1; same <- a$status == b$status; if (same) agree <- agree + 1
    rel <- if (same && a$status == "escaped") abs(a$E/b$E - 1) else NA
    if (!is.na(rel)) worst <- max(worst, rel)
    cat(sprintf("k %3g mu %g |p0| %.3f x0 %.2f: 3D %-15s %-10.6g | brute %-15s %-10.6g %s\n", k, mu,
                sqrt(sum(p0^2)), x0, a$status, a$E, b$status, b$E, if (is.na(rel)) "" else sprintf("rel %.1e", rel)))
  }
}
cat(sprintf("\n%d of %d agree on the outcome; largest relative energy difference %.1e\n", agree, N, worst))
stopifnot(agree == N, worst < 2e-3)   # the brute force has O(dt) error of about 1e-3 at its default step; p4c_converge.R shows it converging onto the tracker
