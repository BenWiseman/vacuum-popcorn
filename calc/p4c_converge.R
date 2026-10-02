# p4c_converge.R: the two cases where p4b's brute force differed most from the 3D tracker, rerun with
# smaller time steps; the brute force should converge onto the tracker.
source("pt3_lib.R")
src <- readLines("p4b_brute3d.R"); i0 <- grep("^brute3 <- function", src); i1 <- grep("^set.seed", src) - 1
eval(parse(text = src[i0:i1]))
set.seed(29); cases <- list()
for (k in c(0.3, 3)) for (mu in c(2, 4)) for (i in 1:6) {
  x0 <- runif(1)^(1/3)*0.97; d <- rnorm(3); r0 <- x0*d/sqrt(sum(d^2)); q <- rnorm(3); p0 <- q/sqrt(sum(q^2))*10^runif(1, -1, 0.5)
  cases[[length(cases) + 1]] <- list(k = k, mu = mu, r0 = r0, p0 = p0) }
last <- c()
for (j in c(11, 23)) { cs <- cases[[j]]; tab <- make_lagtab(cs$k)
  a <- track3(cs$r0, cs$p0, cs$mu, tab, xmin = 0.02)
  cat(sprintf("case %d (k %g mu %g |p0| %.3f): tracker %.7g\n", j, cs$k, cs$mu, sqrt(sum(cs$p0^2)), a$E))
  for (dt in c(5e-5, 1e-5, 2e-6)) { b <- brute3(cs$r0, cs$p0, cs$mu, cs$k, xmin = 0.02, dt0 = dt)
    cat(sprintf("   brute dt %.0e: %.7g  (rel to tracker %.1e)\n", dt, b$E, abs(b$E/a$E - 1))); flush.console() }
  last <- c(last, abs(b$E/a$E - 1)) }
stopifnot(all(last < 1e-4))
