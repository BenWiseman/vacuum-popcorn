# p2_validate.R: the event tracker against brute-force time stepping, where double precision is safe
# (mu <= 5 and a pocket that ends at 1e-2 of its starting radius, so no encounter needs gamma above ~1e5).
source("pt_lib.R")
brute <- function(x0, mu, k, xmin = 1e-2, dt0 = 1e-4) {
  # positions on a line through the centre: particle z (signed), wall at -X and +X
  X <- 1 - 1e-12; t <- 0; z <- x0; E <- 1; pz <- 0; m <- 1   # a wall exactly at rest would never start moving in this stepping
  wallv <- function(X) beta_x(X, k)
  repeat {
    dt <- min(dt0, 0.002*X)
    vz <- pz/E
    # RK2 for the wall
    b1 <- wallv(X); Xm <- X - 0.5*dt*b1; Xn <- X - dt*wallv(Xm)
    zn <- z + vz*dt
    if (Xn < xmin) return(list(status = "trapped_at_end", E = E))
    if (abs(zn) >= Xn) {                         # crossed: bisect for the meeting time
      lo <- 0; hi <- dt
      for (it in 1:60) { mid <- (lo + hi)/2
        Xmid <- X - mid*wallv(X - 0.5*mid*wallv(X)); zmid <- z + vz*mid
        if (abs(zmid) >= Xmid) hi <- mid else lo <- mid }
      Xh <- X - hi*wallv(X - 0.5*hi*wallv(X)); s <- sign(z + vz*hi); if (s == 0) s <- 1
      pn <- s*pz                                 # outward normal is along s
      r <- encounter(E, pn, 0, m, mu, gam(Xh, k))
      if (r$type == "transmit") return(list(status = "escaped", E = r$E))
      if (r$type == "reflect") { E <- r$E; pz <- s*r$pn }
      t <- t + hi; X <- Xh; z <- s*Xh*(1 - 1e-12)
      next
    }
    t <- t + dt; X <- Xn; z <- zn
  }
}
set.seed(5); worst <- 0; agree <- 0; N <- 0
for (mu in c(1.5, 3, 5)) for (k in c(0.3, 3)) for (x0 in c(0.1, 0.3, 0.5, 0.7, 0.9, 0.98)) {
  a <- track(x0, mu, k, xmin = 1e-2); b <- brute(x0, mu, k, xmin = 1e-2)
  N <- N + 1
  same <- a$status == b$status
  rel <- if (same && a$status == "escaped") abs(a$E/b$E - 1) else NA
  if (same) agree <- agree + 1
  if (!is.na(rel)) worst <- max(worst, rel)
  cat(sprintf("mu %4g k %3g x0 %.2f: tracker %-15s E %-10.6g | brute %-15s E %-10.6g %s\n", mu, k, x0,
              a$status, a$E, b$status, b$E, if (is.na(rel)) "" else sprintf("rel %.1e", rel))); flush.console()
}
cat(sprintf("\n%d of %d cases agree on the outcome; largest relative energy difference %.1e\n", agree, N, worst))
stopifnot(agree == N, worst < 1e-3)

# planted failure: a tracker that never lets the wall catch a particle from behind (it always sends the
# particle through the centre) must disagree with the brute force somewhere, or this check proves nothing
track_bad <- function(x0, mu, k, xmin) {
  body <- deparse(track)
  i <- grep("if \\(F0 > 0\\)", body); body[i] <- sub("F0 > 0", "TRUE", body[i])
  f <- eval(parse(text = body)); f(x0, mu, k, xmin = xmin)
}
nbad <- 0
for (mu in c(1.5, 3, 5)) for (k in c(0.3, 3)) for (x0 in c(0.7, 0.9, 0.98)) {
  a <- tryCatch(track_bad(x0, mu, k, 1e-2), error = function(e) list(status = "error", E = NA))
  b <- brute(x0, mu, k, xmin = 1e-2)
  if (a$status != b$status || (a$status == "escaped" && abs(a$E/b$E - 1) > 1e-3)) nbad <- nbad + 1
}
cat(sprintf("planted failure (no catch-ups): %d of 18 cases now disagree with the brute force\n", nbad))
stopifnot(nbad > 0)
