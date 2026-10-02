# p15_schematic.R: the letter's first figure, a sketch of the sequence from pocket to air showers.
# Nothing here is computed; the picture follows the text: a pocket of false vacuum holding trapped
# particles (section 1), its collapse throwing them out (sections 2 to 4), their decay, an unobstructed
# path through the halo, and showers arriving together from one direction (section 5).
blue <- "#1b6ca8"; orange <- "#c2410c"; ink <- "#1f2933"; soft <- "#6b7280"
circ <- function(x, y, r, n = 200) { t <- seq(0, 2*pi, length.out = n); list(x = x + r*cos(t), y = y + r*sin(t)) }
ring <- function(x, y, r, ...) { c1 <- circ(x, y, r); polygon(c1$x, c1$y, ...) }
arr <- function(x0, y0, x1, y1, col = ink, lwd = 1.2, len = 0.05, ...) arrows(x0, y0, x1, y1, length = len, col = col, lwd = lwd, ...)
set.seed(7)
for (dev in c("pdf", "png")) {
  f <- paste0("../letter/fig_schematic.", dev)
  if (dev == "pdf") pdf(f, width = 7.2, height = 2.05) else png(f, width = 1440, height = 410, res = 200)
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new(); plot.window(xlim = c(0, 100), ylim = c(0, 28.5), asp = 1)
  yc <- 14.5; top <- 26.3; bot <- 3.1; ct <- 0.7; cs <- 0.6
  head <- function(x, s) text(x, top, s, font = 2, cex = ct, col = ink)
  foot <- function(x, s) text(x, bot, s, cex = cs, col = soft)

  # 1. the pocket
  x1 <- 9.5; R <- 6.6
  ring(x1, yc, R, col = adjustcolor(blue, 0.14), border = blue, lwd = 2.2)
  i <- 1:19; a <- i*2.39996 + runif(19, -0.25, 0.25); rr <- R*0.84*sqrt((i - 0.5)/19)   # an even scatter
  points(x1 + rr*cos(a), yc + rr*sin(a), pch = 16, cex = 0.42, col = orange)
  head(x1, "Pocket of false vacuum"); foot(x1, "trapped particles:\nlight inside, heavy outside")

  # 2. the pop
  x2 <- 30; r2 <- 2.6
  ring(x2, yc, R, col = NA, border = soft, lwd = 0.8, lty = 3)
  ring(x2, yc, r2, col = adjustcolor(blue, 0.14), border = blue, lwd = 2.2)
  for (th in c(0, 90, 180, 270)*pi/180) arr(x2 + (R - 0.4)*cos(th), yc + (R - 0.4)*sin(th), x2 + (r2 + 0.9)*cos(th), yc + (r2 + 0.9)*sin(th), col = blue, lwd = 1.3, len = 0.045)
  for (th in c(35, 60, 125, 150, 215, 240, 305, 330)*pi/180) arr(x2 + 0.9*cos(th), yc + 0.9*sin(th), x2 + (R + 1.9)*cos(th), yc + (R + 1.9)*sin(th), col = orange, lwd = 1.1, len = 0.04)
  head(x2, "It collapses: a pop")
  text(x2, bot + 0.95, "the closing wall throws them out", cex = cs, col = soft)
  text(x2, bot - 1.05, expression("with energy near 0.01" ~ m[out]^2 / m[`in`]), cex = cs, col = soft)

  # 3. decay
  x3 <- 49
  arr(x3 - 4.4, yc, x3 - 0.6, yc, col = orange, lwd = 1.4)
  points(x3, yc, pch = 16, cex = 0.55, col = ink)
  w <- seq(0, 1, length.out = 120)                               # photon: a wavy line
  lines(x3 + 6.2*w, yc + 4.4*w + 0.45*sin(2*pi*7*w), col = orange, lwd = 1.1)
  segments(x3, yc, x3 + 6.8, yc, col = orange, lwd = 1.1, lty = 2) # neutrino: dashed
  polygon(c(x3, x3 + 6.3, x3 + 5.1), c(yc, yc - 3.5, yc - 5.3), col = adjustcolor(orange, 0.25), border = orange, lwd = 0.9)  # jet
  text(x3 + 7.0, yc + 4.9, expression(gamma), cex = cs + 0.08, col = ink, adj = 0)
  text(x3 + 7.4, yc, expression(nu), cex = cs + 0.08, col = ink, adj = 0)
  text(x3 + 6.7, yc - 4.6, "jet", cex = cs, col = ink, adj = 0)
  head(x3 + 1, "Escapees decay"); foot(x3 + 1, "to photons, neutrinos\nand jets of hadrons")

  # 4. the path
  x4 <- 70.5
  for (dy in c(-1.7, 0, 1.7)) arr(x4 - 7, yc + dy, x4 + 7, yc + dy, col = orange, lwd = 1.1)
  text(x4, yc + 4.6, "the Galaxy's halo", cex = cs, col = soft, font = 3)
  head(x4, "A clear path"); foot(x4, "nothing absorbs them\non the way")

  # 5. Earth
  x5 <- 91.5; g <- 8.6
  segments(x5 - 7.5, g, x5 + 7.5, g, col = ink, lwd = 1.3)
  points(seq(x5 - 6.5, x5 + 6.5, length.out = 9), rep(g, 9), pch = 15, cex = 0.36, col = ink)
  for (xs in c(-3.6, 0.4, 4.4)) {
    xt <- x5 + xs - 5.2; yt <- g + 12.2; xm <- x5 + xs - 1.8; ym <- g + 4.2
    segments(xt, yt, xm, ym, col = orange, lwd = 1.1)
    polygon(c(xm, x5 + xs - 1.5, x5 + xs + 1.5), c(ym, g, g), col = adjustcolor(orange, 0.25), border = orange, lwd = 0.9)
  }
  head(x5, "Earth"); foot(x5, "several air showers at one\ninstant, from one direction")

  # the links between stages
  for (xa in c(18.8, 41.0, 60.6, 80.2)) arr(xa - 1.2, yc, xa + 1.2, yc, col = ink, lwd = 1.6, len = 0.06)
  dev.off()
}
cat("schematic written: five stages, pocket, pop, decay, path, Earth\n")
