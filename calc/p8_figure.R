# p8_figure.R: the letter's figure. (a) escape spectra for a gas at rest, m_out/m_in = 1e6, three wall
# shapes; (b) a hot thermal gas, k = 1, m_out/E_typ = 1e2 and 1e4. Reads p3_spectrum.rds; runs (b).
source("pt3_lib.R")
r3 <- readRDS("p3_spectrum.rds")
pick <- function(mu, k, xmin) Filter(function(r) r$mu == mu && r$k == k && r$xmin == xmin, r3)[[1]]
set.seed(43); tab <- make_lagtab(1)
hot <- lapply(c(1e2, 1e4), function(ratio) {
  Et <- 1000; mu <- ratio*Et
  E <- vapply(1:1000, function(i) {
    x0 <- runif(1)^(1/3)*(1 - 1e-9); d <- rnorm(3); r0 <- x0*d/sqrt(sum(d^2))
    pm <- rgamma(1, shape = 3, scale = Et/3); w <- rnorm(3); p0 <- pm*w/sqrt(sum(w^2))
    r <- track3(r0, p0, mu, tab, xmin = 1e-12); if (r$status == "escaped") r$E else NA }, 0)
  E[!is.na(E)]/(mu^2/Et) })
dens <- function(s) { d <- density(log10(s), bw = 0.12); d }
for (dev in c("pdf", "png")) {
  f <- paste0("../letter/fig_spectra.", dev)
  if (dev == "pdf") pdf(f, width = 7.2, height = 3.1) else png(f, width = 1440, height = 620, res = 200)
  par(mfrow = c(1, 2), mar = c(4, 4, 1.5, 0.6), mgp = c(2.4, 0.7, 0), cex = 0.8)
  cols <- c("#1b6ca8", "#c2410c", "#4d7c0f")
  ds <- lapply(c(0.1, 1, 10), function(k) dens(pick(1e6, k, 1e-12)$s_all))
  plot(NA, xlim = c(-5.5, 1), ylim = c(0, max(sapply(ds, function(d) max(d$y)))*1.05),
       xlab = expression(log[10]~E[out]/(m[out]^2/m[`in`])), ylab = "particles per unit log E", main = "(a) gas at rest")
  ltys <- c(1, 2, 4)      # line types too: the green and the orange are hard to tell apart with red-green colour blindness
  for (i in 1:3) lines(ds[[i]], col = cols[i], lwd = 2, lty = ltys[i])
  legend("topleft", c("k = 0.1", "k = 1", "k = 10"), col = cols, lwd = 2, lty = ltys, bty = "n", seg.len = 3)
  dh <- lapply(hot, dens)
  plot(NA, xlim = c(-5, 2), ylim = c(0, max(sapply(dh, function(d) max(d$y)))*1.05),
       xlab = expression(log[10]~E[out]/(m[out]^2/E[typ])), ylab = "particles per unit log E", main = "(b) hot gas, k = 1")
  for (i in 1:2) lines(dh[[i]], col = cols[i], lwd = 2, lty = ltys[i])
  legend("topleft", c(expression(m[out]/E[typ]==10^2), expression(m[out]/E[typ]==10^4)), col = cols[1:2], lwd = 2, lty = ltys[1:2], bty = "n", seg.len = 3)
  dev.off()
}
cat(sprintf("figure written; hot-gas medians %.3g and %.3g of m_out^2/E_typ\n", median(hot[[1]]), median(hot[[2]])))
