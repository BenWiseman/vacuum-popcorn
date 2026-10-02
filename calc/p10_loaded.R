# p10_loaded.R: a collapsing pocket whose trapped particles slow its wall.
#
# Everything in p3 to p8 assumes a light load: the wall follows gamma(x) = (k (1 - x^3) + 1)/x^2 whatever
# the particles take from it. Here the energy the wall gives is taken off the wall (pl_lib.R, checked in
# p9_loaded_validate.R). The load is
#   load = N (m_out^2 / m_in) / E_pocket,
# the number of trapped particles times m_out^2/m_in, in units of the pocket's own energy (wall plus
# vacuum). A light-load escapee carries about one per cent of m_out^2/m_in, so a load of 100 asks the
# pocket for about as much energy as it has.
#
# Gas at rest, uniform in volume. Energies are quoted as shares of m_out^2/m_in.
#   escaped   share of the particles that leave the pocket
#   stalled   share still inside when the wall has given up all its energy and stops
#   inside    share still inside when the pocket has shrunk to its wall thickness (x_min) with the wall still moving
#   E50, E90  median and 90th percentile of the escapees' energies
#   to esc.   share of the pocket's energy that the escapees carry off
#   to all    share of the pocket's energy given to particles, escaped or not
#   stall at  the pocket's size, as a fraction of its starting size, where the wall stops (0: it never does)
source("pl_lib.R")
args <- commandArgs(trailingOnly = TRUE)
cores <- if (length(args)) as.integer(args[1]) else 12
say <- function(...) { cat(...); if (nzchar(Sys.getenv("P10_PROGRESS"))) cat(..., file = stderr()) }
row <- function(a) say(sprintf("%-6.0e %-4g %-6.0e %-5g | %-7.3f %-7.3f %-7.3f | %-8.4f %-8.4f | %-7.3f %-7.3f | %-8.2g | %2d %-7.1e %-7.1e %6d %4d\n",
  a$mu, a$k, a$xmin, a$ell, a$f_esc, a$f_stalled, a$f_trap, a$q50, a$q90, a$eff, a$drained, a$x_stall,
  a$iters, a$resid, a$overshoot, a$Nsim, a$per_decade))
head <- function() say(sprintf("%-6s %-4s %-6s %-5s | %-7s %-7s %-7s | %-8s %-8s | %-7s %-7s | %-8s | %s\n", "mu", "k", "xmin", "load",
  "escaped", "stalled", "inside", "E50", "E90", "to esc.", "to all", "stall at", "its resid overshoot shells cells/decade"))
out <- list()
solve <- function(mu, k, xmin, ell, per_decade = 60, ...) {
  a <- solve_loaded(mu, k, xmin, ell, cores = cores, per_decade = per_decade, ...)
  a$W <- NULL; a$lam <- NULL; a$per_decade <- per_decade; out[[length(out) + 1]] <<- a
  a
}
loads <- c(0, 3, 10, 20, 40, 80)
say("A. Gas at rest, pocket ending at 1e-6 of its starting size\n"); head()
for (mu in c(1e4, 1e6)) for (k in c(0.1, 1, 10)) { for (ell in loads) row(solve(mu, k, 1e-6, ell)); say("\n") }

say("B. The same for a pocket ending at 1e-12 of its starting size (mu = 1e4, k = 1)\n"); head()
for (ell in loads) row(solve(1e4, 1, 1e-12, ell, s_lo = 1e-7))

# How much do the answers depend on the sampling? Twice the shells, twice the near-wall shells, cells a
# quarter as wide, a quarter of the splitting threshold.
say("\nC. Sampling: mu = 1e4, k = 1, x_min = 1e-6, standard run against a finer one\n"); head()
worst <- c(esc = 0, E50 = 0, eff = 0)
for (ell in c(10, 20, 40, 80)) {
  a <- solve(1e4, 1, 1e-6, ell); row(a)
  b <- solve(1e4, 1, 1e-6, ell, per_decade = 240, N0 = 40000, edge = 100, delta = 5e-4); row(b)
  worst <- pmax(worst, c(abs(a$f_esc - b$f_esc), abs(a$q50/b$q50 - 1), abs(a$eff - b$eff)))
}
say(sprintf("largest change on refining: escaped share %.3f, median energy %.1f per cent, escapees' share of the pocket's energy %.3f\n",
            worst["esc"], 100*worst["E50"], worst["eff"]))
# the checks every row must pass: the wall agrees with its own drain, and no stalled wall has given
# more than one per cent beyond what it had (the fixed point's own mismatch can leave a little)
resid <- vapply(out, `[[`, 0, "resid"); over <- vapply(out, `[[`, 0, "overshoot")
other <- vapply(out, `[[`, 0, "f_other")       # shells whose tracking was cut short: there must be none
say(sprintf("\n%d collapses solved; largest mismatch between a wall and its own drain %.1e of the pocket's energy; largest overshoot %.1e; shells cut short %g\n",
            length(out), max(resid), max(over), max(other)))
saveRDS(lapply(out, function(a) a[c("per_decade", "f_other", "mu", "k", "xmin", "ell", "f_esc", "f_stalled", "f_trap", "q10", "q50", "q90", "q99",
                                    "eff", "drained", "x_stall", "iters", "resid", "overshoot", "Nsim", "f_tail", "e_tail")]),
        "p10_loaded.rds")
stopifnot(max(resid) < 5e-3, max(over) < 0.01, max(other) == 0)
