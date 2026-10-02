# p9_loaded_validate.R: the loaded-collapse solver (pl_lib.R) against three independent checks.
#  1. With no load it must reproduce the light-load tracker (pt_lib.R) shell for shell.
#  2. With a load it must agree with brute-force time stepping, in which every shell and the wall
#     advance together and each encounter drains the wall at once (no grid, no iteration).
#  3. Planted failures: a solver with the wrong load, and one that never drains the wall, must
#     both disagree with the brute force, or check 2 proves nothing.
source("pt_lib.R")
track_light <- track; gam_light <- gam
source("pl_lib.R")

# ---- 1. no load
set.seed(17); n <- 0; agree <- 0; worst <- 0
for (cfg in list(c(1e2, 0.1, 1e-6), c(1e4, 1, 1e-6), c(1e4, 1, 1e-12), c(1e6, 10, 1e-12), c(3, 0.3, 1e-2))) {
  mu <- cfg[1]; k <- cfg[2]; xmin <- cfg[3]
  g <- make_grid(xmin); W <- make_wall(k, g, numeric(length(g)))
  for (x in pmin(runif(300)^(1/3), 1 - 1e-9)) {
    a <- track_light(x, mu, k, xmin = xmin, max_events = 2000); b <- track_loaded(x, mu, W, xmin = xmin)
    n <- n + 1
    if (a$status == b$status) { agree <- agree + 1
      if (a$status == "escaped") worst <- max(worst, abs(a$E/b$E - 1)) }
  }
}
cat(sprintf("1. no load: %d of %d shells have the same fate as the light-load tracker; largest escape-energy difference %.1e\n",
            agree, n, worst))
stopifnot(agree == n, worst < 1e-8)

# ---- 2. brute force with a load
brute_loaded <- function(x0, mu, k, ell, xmin = 1e-2, dt0 = 1e-4) {
  Ns <- length(x0); scale <- ell*(1 + k)/(Ns*mu^2)
  X <- 1 - 1e-12; lam <- 0                    # a wall exactly at rest never starts under this stepping
  z <- x0; E <- rep(1, Ns); pz <- rep(0, Ns); st <- rep("inside", Ns)
  num <- function(X) { o <- 1 - X; o*(1 + X) + k*o*(1 + X + X^2) - lam }
  bet <- function(X) { h <- num(X)/X^2; if (h <= 0) 0 else sqrt(h*(h + 2))/(h + 1) }
  step <- function(X, dt) X - dt*bet(X - 0.5*dt*bet(X))
  finish <- function(why) { st[st == "inside"] <<- why; list(status = st, E = E, lam = lam, X = X) }
  repeat {
    if (lam > 0 && num(X) <= 0) return(finish("stalled"))
    ins <- which(st == "inside"); if (!length(ins)) return(finish("none"))
    dt <- min(dt0, 0.002*X); v <- pz/E
    Xn <- step(X, dt)
    if (Xn < xmin) return(finish("trapped_at_end"))
    hit <- ins[abs(z[ins] + v[ins]*dt) >= Xn]
    if (!length(hit)) { X <- Xn; z[ins] <- z[ins] + v[ins]*dt; next }
    th <- vapply(hit, function(i) { lo <- 0; hi <- dt
      for (it in 1:60) { mid <- (lo + hi)/2
        if (abs(z[i] + v[i]*mid) >= step(X, mid)) hi <- mid else lo <- mid }
      hi }, 0)
    i <- hit[which.min(th)]; t <- min(th)
    Xh <- step(X, t); z[ins] <- z[ins] + v[ins]*t; X <- Xh
    s <- sign(z[i]); if (s == 0) s <- 1
    r <- encounter(E[i], s*pz[i], 0, 1, mu, num(Xh)/Xh^2 + 1)
    if (r$type == "miss") { z[i] <- s*Xh*(1 - 1e-12); next }
    lam <- lam + scale*(r$E - E[i])
    E[i] <- r$E
    if (r$type == "transmit") st[i] <- "escaped" else { pz[i] <- s*r$pn; z[i] <- s*Xh*(1 - 1e-12) }
  }
}
compare <- function(mu, k, ell, x0, xmin = 1e-2, solver = solve_events, dt0 = 1e-4) {
  b <- brute_loaded(x0, mu, k, ell, xmin = xmin, dt0 = dt0)
  a <- solver(mu, k, xmin, ell, x0)
  same <- a$status == b$status; both <- same & a$status == "escaped"
  list(n = length(x0), same = sum(same), worst = if (any(both)) max(abs(a$E[both]/b$E[both] - 1)) else 0,
       drained_a = a$drained, drained_b = b$lam/(1 + k), n_esc = sum(b$status == "escaped"),
       stalled = any(b$status == "stalled"))
}
x0 <- shells(24)
# Loads chosen so that shells escape while the wall is drained, some with a stall and some without:
# a first set of six cases all stalled before any shell escaped, so every fate agreed trivially.
cases <- list(c(1.5, 0.3, 0.3), c(1.5, 3, 1), c(3, 0.3, 0.3), c(3, 3, 0.3), c(5, 3, 0.1), c(5, 0.3, 1),
              c(10, 3, 0.3), c(3, 0.3, 1), c(5, 3, 0.3), c(10, 3, 1))
tot <- 0; ok <- 0; worst <- 0; res <- list()
for (cs in cases) {
  r <- compare(cs[1], cs[2], cs[3], x0); res[[length(res) + 1]] <- r
  tot <- tot + r$n; ok <- ok + r$same; worst <- max(worst, r$worst)
  cat(sprintf("   mu %3g k %3g load %4g: %2d of %d shells agree, %2d escape, energies to %.1e; share of pocket energy given to particles %.4f (solver) %.4f (brute)%s\n",
              cs[1], cs[2], cs[3], r$same, r$n, r$n_esc, r$worst, r$drained_a, r$drained_b, if (r$stalled) "; wall stalls" else ""))
  flush.console()
}
cat(sprintf("2. with a load: %d of %d shells have the same fate as the brute force; largest escape-energy difference %.1e\n", ok, tot, worst))
dmax <- max(vapply(res, function(r) abs(r$drained_a - r$drained_b), 0))
cat(sprintf("   largest difference in the share of pocket energy given to particles: %.1e\n", dmax))
nesc <- sum(vapply(res, `[[`, 0, "n_esc")); nstall <- sum(vapply(res, `[[`, TRUE, "stalled"))
cat(sprintf("   %d escapes compared in all; the wall stalls in %d of %d cases\n", nesc, nstall, length(cases)))
stopifnot(ok == tot, worst < 2e-3, dmax < 2e-3, nesc >= 50, nstall >= 2, nstall <= length(cases) - 4)

# ---- 3. planted failures
half <- function(mu, k, xmin, ell, x0) solve_events(mu, k, xmin, ell/2, x0)
none <- function(mu, k, xmin, ell, x0) solve_events(mu, k, xmin, 0, x0)
bad_half <- 0; bad_none <- 0
for (cs in cases) {
  for (nm in c("half", "none")) {
    r <- compare(cs[1], cs[2], cs[3], x0, solver = get(nm))
    wrong <- r$same < r$n || r$worst > 2e-3
    if (nm == "half") bad_half <- bad_half + wrong else bad_none <- bad_none + wrong
  }
}
cat(sprintf("3. planted failures: a solver given half the load disagrees with the brute force in %d of %d cases; one that never drains the wall, in %d of %d\n",
            bad_half, length(cases), bad_none, length(cases)))
stopifnot(bad_half > 0, bad_none > 0)

# ---- 4. the fast solver against the exact one
# solve_events is exact and slow (its work grows as the square of the number of shells); p10 uses
# solve_loaded, which iterates a wall history on a grid. Give both the same 300 equal shells, at a
# mass ratio and loads like p10's, and they must agree shell for shell. The energies compared are those
# of every shell, escaped or not: an escapee's energy hardly depends on how fast the wall was (which is
# the physics of p10), so on escapees alone a solver that never drained its wall would pass.
mu <- 1e4; k <- 1; xmin <- 1e-6; N <- 300; B <- list(a = (seq_len(N) - 1)/N, b = seq_len(N)/N)
same <- 0; tot <- 0; worst <- 0; dmax <- 0
for (ell in c(2, 5, 10)) {
  e <- solve_events(mu, k, xmin, ell, shells(N))
  a <- solve_loaded(mu, k, xmin, ell, base = B, delta = Inf, tol = 1e-4, maxit = 100)
  same <- same + sum(a$status == e$status); tot <- tot + N
  worst <- max(worst, abs(a$E/e$E - 1)); dmax <- max(dmax, abs(a$drained - e$drained))
  cat(sprintf("   load %2g: %d of %d shells agree, %d escape; share of pocket energy given to particles %.4f (fast) %.4f (exact)\n",
              ell, sum(a$status == e$status), N, sum(e$status == "escaped"), a$drained, e$drained)); flush.console()
}
cat(sprintf("4. fast solver against exact: %d of %d shells have the same fate; largest energy difference, any shell, %.1e; largest difference in energy given %.1e\n",
            same, tot, worst, dmax))
stopifnot(same == tot, worst < 1e-3, dmax < 1e-4)
# planted failure: the fast solver stopped after one pass, before the wall has heard about its drain
a1 <- solve_loaded(mu, k, xmin, 10, base = B, delta = Inf, maxit = 1, tol = 1e-4)
gap <- max(abs(a1$E/e$E - 1))
cat(sprintf("   planted failure: stopped after one pass, the fast solver's energies are off by up to %.2f and the energy given by %.4f\n",
            gap, abs(a1$drained - e$drained)))
stopifnot(gap > 0.05, abs(a1$drained - e$drained) > 1e-3)
