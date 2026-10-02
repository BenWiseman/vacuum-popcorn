# pl_lib.R: the radial tracker of pt_lib.R, with a wall that the trapped particles slow down.
#
# Units as in pt_lib.R: c = 1, R0 = 1, particle energies in m_in, mu = m_out/m_in. Energy in wall +
# vacuum + particles is conserved; in units of 4 pi sigma R0^2 the wall holds gamma x^2 and the vacuum
# k x^3, so
#   gamma(x) = (1 + k - lam(x) - k x^3) / x^2,
# where lam(x) is the energy the wall has given to particles by the time it reaches radius x. The wall
# only moves inward until it stops, so x orders the encounters in time.
#
# The particles are spherical shells (a particle at rest on a sphere of radius x0 stands for every
# particle on that sphere), so the wall stays spherical. A shell of weight w stands for a share w of
# the N particles. The load is
#   ell = N (m_out^2 / m_in) / E_pocket,     E_pocket = 4 pi sigma R0^2 (1 + k),
# the energy N particles would carry at m_out^2/m_in each, in units of the pocket's own. An encounter
# that changes a particle's energy by dE (in m_in) drains ell (1 + k) w dE / mu^2 from the wall.
#
# Two solvers. solve_events is exact and slow: it takes the encounters in order and drains the wall
# at each. solve_loaded is fast: it iterates the wall's history, holding lam on the cells of a fixed
# grid. lam[j] applies while the wall is in cell j = [lo[j], hi[j]), and an encounter in cell j drains
# the wall from cell j - 1 inward, a one-cell delay (p9 checks solve_loaded against solve_events shell
# for shell; p10 repeats its runs on narrower cells). Within a cell gamma falls as x grows, so the
# wall's lag behind light is monotone inside a cell; across cells it need not be, since a drained
# wall can slow down while it shrinks. The catch-up search in next_meet is built on that: it scans the
# cell edges and solves inside the first cell where the wall has caught the shell, which is exact
# whatever the wall does between cells.
source("pk_lib.R")

make_grid <- function(xmin, per_decade = 60) {
  lo_x <- 10^seq(log10(xmin) - 1, log10(0.5), length.out = ceiling(per_decade*(log10(0.5) - log10(xmin) + 1)) + 1)
  hi_s <- 10^seq(log10(0.5), -10, length.out = ceiling(per_decade*(10 + log10(0.5))) + 1)   # s = 1 - x
  sort(unique(c(lo_x, 1 - hi_s, 1)))       # cell tops; the last is 1
}

# integrate() to 1e-12, and where it reports that roundoff stops it short of that (very narrow cells),
# its best estimate, which is still far tighter than anything quoted.
quad <- function(f, a, b) {
  r <- integrate(f, a, b, rel.tol = 1e-12, abs.tol = 0, subdivisions = 2000L, stop.on.error = FALSE)
  if (!(r$message %in% c("OK", "roundoff error was detected", "roundoff error is detected in the extrapolation table")))
    stop("quad: ", r$message)
  r$value
}

# reuse: the wall this one was made from by adding one cell edge inside its bottom cell; every cell
# above that is unchanged, so its integral is carried over instead of being computed again.
make_wall <- function(k, g, lam, reuse = NULL) {
  M <- length(g); stopifnot(g[M] == 1, length(lam) == M, all(diff(g) > 0))
  lo <- c(0, g[-M]); hi <- g
  cell <- function(x) pmin(findInterval(x, c(0, g)), M)
  num <- function(x, j, omx = 1 - x) omx*(1 + x) + k*omx*(1 + x + x^2) - lam[j]   # (gamma - 1) x^2
  gm1 <- function(x, j, omx = 1 - x) num(x, j, omx)/x^2
  lagf <- function(x, j, omx = 1 - x) {
    h <- gm1(x, j, omx); g1 <- h + 1; b <- sqrt(h*(h + 2))/g1; 1/(g1^2*(1 + b)*b) }
  # Where the wall's energy runs out: the first cell, going inward, whose top it cannot reach.
  top_num <- num(hi, seq_len(M), 1 - hi); top_num[M] <- 1      # the wall starts at rest at x = 1
  dead <- which(top_num < 0)
  x_stall <- if (length(dead)) hi[max(dead)] else 0
  # Integral of lag over [a, b] inside cell j, in the variable t with x = hi[j] - t^2. A wall that
  # enters a cell with almost no kinetic energy left has a lag like 1/sqrt(distance below the cell's
  # top), as the wall starting from rest has below x = 1; the substitution removes that end point.
  # (Integrating in x itself stopped a run with "extremely bad integrand behaviour" one cell above
  # a stall.)
  oh <- 1 - hi
  piece <- function(a, b, j) {
    if (b <= a) return(0)
    quad(function(t) lagf(hi[j] - t^2, j, oh[j] + t^2)*2*t, sqrt(hi[j] - b), sqrt(hi[j] - a))
  }
  live <- which(lo >= x_stall)                 # cells the wall actually crosses
  P <- rep(NA_real_, M)
  if (!is.null(reuse) && length(reuse$g) == M - 1 && M > 2 && all(reuse$g[-1] == g[-(1:2)]) &&
      all(reuse$lam[-1] == lam[-(1:2)])) { P[3:M] <- reuse$P[-1]; live <- intersect(live, 1:2) }
  for (j in live) P[j] <- piece(lo[j], hi[j], j)
  I_lag <- function(a, b) {                    # integral of (1/beta - 1) from a to b
    if (b <= a) return(0)
    ja <- cell(a); jb <- cell(b)
    if (lo[jb] == b && jb > 1) jb <- jb - 1      # an upper limit on a cell edge closes the cell below
    if (ja == jb) return(piece(a, b, ja))
    mid <- if (jb - ja > 1) sum(P[(ja + 1):(jb - 1)]) else 0
    piece(a, hi[ja], ja) + mid + piece(lo[jb], b, jb)
  }
  # I_lag(lo[i], x) for cells i = cell(x) down to jend + 1, summed outward-in so nothing cancels
  I_from_edges <- function(x, jend, jx = cell(x)) {
    top <- piece(lo[jx], x, jx)
    if (jx <= jend + 1) return(list(j = jx, I = top))
    list(j = jx:(jend + 1), I = top + c(0, cumsum(P[(jx - 1):(jend + 1)])))
  }
  list(k = k, g = g, lo = lo, hi = hi, lam = lam, P = P, cell = cell, gm1 = gm1, lagf = lagf,
       gam = function(x) gm1(x, cell(x)) + 1, I_lag = I_lag, I_from_edges = I_from_edges,
       x_stall = x_stall)
}

# A shell reflected at radius x now moves inward with 1/v - 1 = delta. Where does it next meet the
# wall, above radius xend? Returns list(x, dir): dir "inward" if the wall catches it from behind,
# "out" if it crosses the centre and meets the far side head-on; NULL if neither happens above xend.
next_meet <- function(x, delta, W, xend) {
  if (x <= xend) return(NULL)
  # Ff(X) = I_lag(X, x) - (x - X) delta is the wall's lead over the shell in reaching radius X; the
  # first X < x where it is <= 0 is the catch-up. Ff starts at 0 and grows, because the reflected
  # shell outruns the wall.
  Ff <- function(X) W$I_lag(X, x) - (x - X)*delta
  jx <- W$cell(x); jend <- W$cell(xend); xc <- NA
  if (W$lo[jx] == x && jx > 1) jx <- jx - 1       # reflected exactly on a cell's lower edge: it flies in the cell below
  if (!(W$lagf(x, jx) > delta)) stop("reflected shell is not outrunning the wall")
  bot <- max(W$lo[jx], xend)
  if (W$lagf(bot, jx) < delta) {                  # the wall overtakes the shell's speed inside this cell
    xp <- uniroot(function(X) W$lagf(X, jx) - delta, c(bot, x), tol = 1e-15)$root
    if (Ff(bot) <= 0) xc <- uniroot(Ff, c(bot, xp), tol = 1e-15)$root
  }
  if (is.na(xc) && jx > jend) {
    ed <- W$I_from_edges(x, jend, jx)             # edges lo[jx], ..., lo[jend + 1]
    Xe <- W$lo[ed$j]; Fe <- ed$I - (x - Xe)*delta
    hit <- which(Fe[-1] <= 0)                     # first edge below lo[jx] where the wall has caught up
    if (length(hit)) {
      i <- ed$j[hit[1] + 1]                       # caught inside cell i
      xc <- uniroot(Ff, c(W$lo[i], W$hi[i]), tol = 1e-15)$root
    } else if (Ff(xend) <= 0) {                   # caught in the last cell, above xend
      xc <- uniroot(Ff, c(xend, W$hi[jend]), tol = 1e-15)$root
    }
  }
  if (!is.na(xc)) return(list(x = xc, dir = "inward"))
  # Not caught above xend: does it cross the centre and meet the far side of the wall above xend?
  # H(rho) = I_lag(rho, x) - 2 rho - (x + rho) delta is zero where they meet, and falls as rho grows.
  H <- function(rho) W$I_lag(rho, x) - 2*rho - (x + rho)*delta
  if (H(xend) <= 0) return(NULL)                  # still inside when the wall reaches xend
  up <- min(x, xend + W$I_lag(xend, x)/2)         # 2 rho <= I_lag(rho, x) <= I_lag(xend, x)
  list(x = uniroot(H, c(xend, up), tol = 1e-14*up)$root, dir = "out")
}

# One shell, initially at rest at x0. Returns its fate and every encounter as (radius, energy change).
track_loaded <- function(x0, mu, W, xmin = 1e-6, max_events = 50000) {
  xend <- max(xmin, W$x_stall)
  end_status <- if (W$x_stall > xmin) "stalled" else "trapped_at_end"
  E <- 1; p <- 0; A <- 1; x <- x0; dir <- "rest"; n <- 0; nc <- 0
  ex <- numeric(0); edE <- numeric(0)
  res <- function(status, E) list(status = status, E = E, n = n, x = x, ncatch = nc, ex = ex, edE = edE)
  if (x0 < xend) return(res(end_status, E))
  repeat {
    n <- n + 1
    if (n > max_events) return(res("too_many", E))
    g <- W$gam(x)
    pn <- switch(dir, rest = 0, out = p, inward = -p)
    if (dir == "inward") nc <- nc + 1
    r <- encounter(E, pn, 0, 1, mu, g, A = A)
    if (r$type == "miss") stop("a shell missed the wall, which the geometry forbids")
    if (r$E > 2*g*mu*(1 + 1e-9)) stop("encounter bound E_after < 2 g m_out violated")
    ex <- c(ex, x); edE <- c(edE, r$E - E)
    if (r$type == "transmit") return(res("escaped", r$E))
    E <- r$E; p <- -r$pn; A <- r$A                  # now moving inward from x; A/p = 1/v - 1
    nx <- next_meet(x, A/p, W, xend)
    if (is.null(nx)) return(res(end_status, E))
    x <- nx$x; dir <- nx$dir
    if (dir == "out") A <- E + p
  }
}

# Shells at the midpoints of N equal volumes: every shell stands for the same number of particles.
shells <- function(N) ((seq_len(N) - 0.5)/N)^(1/3)
wquant <- function(v, w, p) {                 # weighted quantile
  if (!length(v)) return(NA_real_)
  o <- order(v); cw <- cumsum(w[o])/sum(w); v[o][which(cw >= p)[1]]
}

# Fates of weighted shells. w sums to 1; a shell of weight w stands for a share w of the particles.
summarise_fates <- function(st, E, w, mu, ell, k) {
  esc <- st == "escaped"; s <- E[esc]/mu^2; ws <- w[esc]
  list(mu = mu, k = k, ell = ell, Nsim = length(st), status = st, E = E, w = w, s = s, ws = ws,
       f_esc = sum(ws), f_trap = sum(w[st == "trapped_at_end"]), f_stalled = sum(w[st == "stalled"]),
       f_other = sum(w[!(st %in% c("escaped", "trapped_at_end", "stalled"))]),   # shells cut off at max_events
       q10 = wquant(s, ws, 0.1), q50 = wquant(s, ws, 0.5), q90 = wquant(s, ws, 0.9), q99 = wquant(s, ws, 0.99),
       eff = ell*sum(ws*(E[esc] - 1))/mu^2,                 # energy the escapees took / pocket energy
       f_tail = sum(ws[s > 1])/max(sum(ws), 1e-300),        # share of escapees above m_out^2/m_in
       e_tail = sum((ws*s)[s > 1])/max(sum(ws*s), 1e-300))  # their share of the escapees' energy
}

# Fixed-point solver for many light shells: wall history -> tracks -> drain -> wall history, until the
# drain the tracks produce is the drain the wall was given. alpha damps the update. lam lives on the
# cells of a fixed grid, so a drain acts from the cell below the one it happened in.
#
# A few shells take far more than their share: a shell the wall catches from behind at Lorentz factor
# g leaves with up to 2 g m_out, and that energy varies steeply with where the shell started. A shell
# standing for 1/N of the particles then overstates what its neighbours take. So any shell whose
# largest single encounter drains more than delta of the pocket's energy is split into `split`
# thinner shells, again and again, until no encounter drains more than delta. So is any shell whose
# largest encounter differs from a neighbour's by that much: splitting only the shells that look heavy
# would thin out the peaks that were sampled and leave unsampled the peaks that fell between shells,
# and the total would come out low. Shells are intervals [a, b] of the volume fraction u = x0^3, so
# the weight is b - a.
#
# The shells that matter most for the wall's energy start close to it. The wall, starting from rest,
# pushes them ahead of itself and catches them again and again, and a catch-up at Lorentz factor g
# can hand over up to 2 g m_out. Equal-volume shells put almost none there, so the starting set adds
# shells spaced evenly in log(1 - x0), `edge` per decade, from s_lo (the wall's own thickness, or
# 1e-7 if that is smaller) out to 0.1.
base_shells <- function(N0, edge, s_lo) {
  u <- seq(0, 1, length.out = N0 + 1)
  if (edge > 0) u <- c(u, (1 - 10^seq(log10(s_lo), -1, length.out = ceiling(edge*(-1 - log10(s_lo))) + 1))^3)
  u <- sort(unique(u)); list(a = u[-length(u)], b = u[-1])
}
solve_loaded <- function(mu, k, xmin, ell, N0 = 20000, edge = 50, s_lo = max(xmin, 1e-7), per_decade = 60,
                         delta = 2e-3, split = 4, w_min = 1e-10, max_shells = 1.5e5, alpha = 0.5, maxit = 60,
                         tol = 2e-3, cores = 1, base = base_shells(N0, edge, s_lo), verbose = FALSE) {
  ua <- base$a; ub <- base$b
  g <- make_grid(xmin, per_decade); M <- length(g); lam <- numeric(M)
  run_tracks <- function(x0, W) {
    f <- function(x) track_loaded(x, mu, W, xmin = xmin)
    tr <- if (cores > 1 && length(x0) > 50) parallel::mclapply(x0, f, mc.cores = cores) else lapply(x0, f)
    bad <- vapply(tr, function(t) inherits(t, "try-error") || is.null(t), TRUE)
    if (any(bad)) stop(sprintf("%d shells failed; first: %s", sum(bad), as.character(tr[[which(bad)[1]]])))
    tr
  }
  top <- function(tr) vapply(tr, function(t) if (length(t$edE)) max(t$edE) else 0, 0)
  heavy <- function(tr) {                     # shells kept in order of ua, so neighbours are adjacent
    tp <- top(tr); n <- length(tp)
    jump <- pmax(abs(tp - c(tp[1], tp[-n])), abs(tp - c(tp[-1], tp[n])))
    h <- ell*(ub - ua)*pmax(tp, jump)/mu^2
    big <- which(h > delta & (ub - ua) > w_min)
    room <- floor((max_shells - length(ua))/(split - 1))     # out of room: split the heaviest first
    if (length(big) > room) big <- big[order(-h[big])][seq_len(max(room, 0))]
    big
  }
  sweep <- function(W) {                      # tracks for every shell, splitting the heavy ones
    tr <- run_tracks(((ua + ub)/2)^(1/3), W)
    repeat {
      big <- heavy(tr)
      if (!length(big)) break
      a <- rep(ua[big], each = split); h <- rep((ub[big] - ua[big])/split, each = split)
      na <- a + h*(seq_len(split) - 1); nb <- na + h
      a2 <- c(ua[-big], na); o <- order(a2)
      tr <- c(tr[-big], run_tracks(((na + nb)/2)^(1/3), W))[o]
      ua <<- a2[o]; ub <<- c(ub[-big], nb)[o]
    }
    list(tr = tr, left = sum(ell*(ub - ua)*top(tr)/mu^2 > delta))
  }
  drain <- function(tr, W) {                  # lam[j] = (energy given in cells above j), in 4 pi sigma R0^2
    n <- vapply(tr, function(t) length(t$ex), 0L)
    allx <- unlist(lapply(tr, `[[`, "ex")); alld <- unlist(lapply(tr, `[[`, "edE"))*rep(ub - ua, n)
    D <- ell*(1 + k)/mu^2*as.numeric(rowsum(c(alld, numeric(M)), c(W$cell(allx), seq_len(M)))[, 1])
    list(lam = rev(cumsum(rev(D))) - D, cell = D)
  }
  it <- 0; trail <- numeric(0)
  if (ell > 0) for (it in seq_len(maxit)) {
    W <- make_wall(k, g, lam); sw <- sweep(W); d <- drain(sw$tr, W)
    change <- max(abs(d$lam - lam))/(1 + k); trail <- c(trail, change)
    if (verbose) cat(sprintf("  it %2d change %.2e drained %.4f stall %.3g shells %d\n", it, change,
                             d$lam[1]/(1 + k), W$x_stall, length(ua)))
    if (change < tol) break
    lam <- (1 - alpha)*lam + alpha*d$lam
  }
  W <- make_wall(k, g, lam); sw <- sweep(W); tr <- sw$tr; d <- drain(tr, W)
  st <- vapply(tr, `[[`, "", "status"); E <- vapply(tr, `[[`, 0, "E")
  xs <- W$x_stall; jend <- W$cell(max(xmin, xs))
  drained <- (d$lam[jend] + d$cell[jend])/(1 + k)
  if (xs > 0) {
    # The wall stops somewhere inside cell jend, the last one it enters. Every encounter there was
    # worked out with the energy the wall had on entering the cell, so together they can take more
    # than it has. Take them in the order the wall meets them and stop the wall at the first one it
    # has not the energy for; that encounter and the later ones did not happen, so the shells go back to the
    # energy they had and are counted as stalled.
    n <- vapply(tr, function(t) length(t$ex), 0L); sh <- rep(seq_along(tr), n)
    allx <- unlist(lapply(tr, `[[`, "ex")); dE <- unlist(lapply(tr, `[[`, "edE"))
    o <- which(W$cell(allx) == jend); o <- o[order(-allx[o])]
    given <- cumsum(ell*(1 + k)/mu^2*dE[o]*(ub - ua)[sh[o]])
    cut <- which(given > W$gm1(allx[o], jend)*allx[o]^2)[1]   # (gamma - 1) x^2: the wall's kinetic energy
    if (!is.na(cut)) {
      gone <- o[cut:length(o)]; xs <- allx[o[cut]]
      back <- rowsum(dE[gone], sh[gone])
      i <- as.integer(rownames(back)); E[i] <- E[i] - back[, 1]; st[i] <- "stalled"
      drained <- (d$lam[jend] + if (cut > 1) given[cut - 1] else 0)/(1 + k)
    }
  }
  avail <- ((1 - xs^2) + k*(1 - xs^3))/(1 + k)            # what a wall stopped at xs has given up
  c(summarise_fates(st, E, ub - ua, mu, ell, k),
    list(xmin = xmin, iters = it, trail = trail,
         resid = max(abs(d$lam - lam))/(1 + k),             # how far the final wall is from its own drain
         x_stall = xs, lam = lam, W = W, ua = ua, ub = ub,
         drained = drained,                                 # all energy given to particles / pocket energy
         overshoot = if (xs > 0) drained - avail else max(0, drained - 1),
         heaviest = max(ell*(ub - ua)*top(tr)/mu^2),        # largest single encounter, in pocket energies
         unsplit = sw$left,
         worst_cell = max(d$cell)/(1 + k)))                 # largest single-cell drain: the delay's size
}

# Event-ordered solver, exact for any number of shells and slow for many. The wall only moves inward,
# so an encounter at radius x can depend only on encounters at larger radii. Take the pending
# encounter at the largest radius, apply it, drain the wall from that radius inward, and work out
# again where every shell in flight will next meet the slowed wall. No grid and no iteration: lam is
# a step function with a step at each encounter.
solve_events <- function(mu, k, xmin, ell, x0, w = rep(1/length(x0), length(x0)), max_events = 200000,
                         verbose = FALSE) {
  Ns <- length(x0); scale <- ell*(1 + k)*w/mu^2
  E <- rep(1, Ns); p <- rep(0, Ns); Am <- rep(1, Ns)      # Am = E - p, kept without cancellation
  xl <- x0                                    # where each shell last left the wall
  xn <- x0; nd <- rep("rest", Ns)             # where and how it next meets the wall
  st <- rep("inside", Ns); nev <- integer(Ns)
  edges <- 1; lams <- 0                       # cell tops, descending, and the lam of the cell under each
  W <- make_wall(k, 1, 0); drained <- 0; n <- 0
  finish <- function(why) {
    st[st == "inside"] <- why
    c(summarise_fates(st, E, w, mu, ell, k),
      list(xmin = xmin, x_stall = W$x_stall, drained = drained/(1 + k), events = n))
  }
  repeat {
    ins <- which(st == "inside" & !is.na(xn))
    if (!length(ins)) return(finish("trapped_at_end"))
    i <- ins[which.max(xn[ins])]; x <- xn[i]
    if (x < xmin) return(finish("trapped_at_end"))
    n <- n + 1; if (n > max_events) stop("solve_events: too many events")
    g <- W$gm1(x, 1) + 1                      # every earlier encounter was at a larger radius: bottom cell
    r <- if (nd[i] == "out") encounter(E[i], p[i], 0, 1, mu, g, A = E[i] + p[i])
         else encounter(E[i], -p[i], 0, 1, mu, g, A = Am[i])
    if (r$type == "miss") stop("a shell missed the wall, which the geometry forbids")
    drained <- drained + scale[i]*(r$E - E[i]); E[i] <- r$E; nev[i] <- nev[i] + 1L
    if (r$type == "transmit") { st[i] <- "escaped"; xn[i] <- NA }
    else { p[i] <- -r$pn; Am[i] <- r$A; xl[i] <- x; nd[i] <- "flying" }
    if (x < edges[length(edges)]) { edges <- c(edges, x); lams <- c(lams, drained) }
    else lams[length(lams)] <- drained        # two encounters at one radius
    W <- make_wall(k, rev(edges), rev(lams), reuse = W)
    if (verbose) cat(sprintf("  event %4d shell %3d at x %.6g: %s, drained %.4f\n", n, i, x, r$type, drained/(1 + k)))
    if (W$x_stall > 0) return(finish("stalled"))
    for (j in which(st == "inside" & nd != "rest")) {    # every shell in flight now meets a slower wall
      nx <- next_meet(xl[j], Am[j]/p[j], W, xmin)
      if (is.null(nx)) { xn[j] <- NA; nd[j] <- "flying" } else { xn[j] <- nx$x; nd[j] <- nx$dir }
    }
  }
}
