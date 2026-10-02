# p1_kinematics.R: tests of pk_lib.R against closed forms and against 100-digit arithmetic.
source("pk_lib.R")
set.seed(3)
cat("1. CLOSED FORMS\n")
g <- 37; b <- sqrt(1 - 1/g^2)
r <- encounter(1, 0, 0, 1, 1e6, g)          # at rest, reflected
cat(sprintf("   at rest, reflected: E/m = %.10g, closed form g^2 (1 + b^2) = %.10g\n", r$E, g^2*(1 + b^2)))
stopifnot(r$type == "reflect", abs(r$E/(g^2*(1 + b^2)) - 1) < 1e-13)
r <- encounter(5, 5, 0, 0, 1e9, g)          # massless, head-on
cat(sprintf("   massless head-on, reflected: E_after/E = %.10g, closed form g^2 (1 + b)^2 = %.10g\n", r$E/5, g^2*(1 + b)^2))
stopifnot(abs(r$E/5/(g^2*(1 + b)^2) - 1) < 1e-13)
# the product rule for reflection: E_before E_after = g^2 (m_in^2 + pt2) + u^2
for (i in 1:2000) {
  m_in <- 10^runif(1, -3, 0); m_out <- m_in*10^runif(1, 0.01, 6); pt2 <- (m_in*10^runif(1, -3, 1))^2
  E <- sqrt(m_in^2 + pt2 + (m_in*10^runif(1, -3, 4))^2)
  pn <- sample(c(-1, 1), 1)*sqrt(E^2 - m_in^2 - pt2); gg <- 1 + 10^runif(1, -4, 4)
  r <- encounter(E, pn, pt2, m_in, m_out, gg)
  if (r$type == "reflect") stopifnot(abs(E*r$E/(gg^2*(m_in^2 + pt2) + r$u^2) - 1) < 1e-9)
}
cat("   product rule E_before E_after = g^2 (m_in^2 + pt^2) + u^2: holds on 2000 random reflections\n")
cat("\n2. THE SWEEP LIMIT: a particle at rest, wall Lorentz factor rising (m_out/m_in = 1e6)\n")
for (gg in c(1e7, 1e9, 1e12, 1e15)) {
  r <- encounter(1, 0, 0, 1, 1e6, gg)
  cat(sprintf("   g = %.0e: %s, E/m_in = %.12g, limit (m_out^2 + m_in^2)/(2 m_in) = %.12g\n",
              gg, r$type, r$E, sweep_limit(1, 0, 1e6)))
}
stopifnot(abs(encounter(1, 0, 0, 1, 1e6, 1e15)$E/sweep_limit(1, 0, 1e6) - 1) < 1e-9)
# write random cases for the 100-digit reference, covering huge g and huge energies
n <- 3000
cases <- data.frame(m_in = 10^runif(n, -6, 0))
cases$m_out <- cases$m_in*10^runif(n, 0.01, 9)
cases$pt2 <- (cases$m_in*10^runif(n, -4, 2))^2*(runif(n) < 0.7)
cases$p <- cases$m_in*10^runif(n, -4, 12)
cases$sgn <- sample(c(-1, 1), n, replace = TRUE)
cases$g <- 1 + 10^runif(n, -6, 14)
cases$E <- sqrt(cases$m_in^2 + cases$pt2 + cases$p^2)
cases$pn <- cases$sgn*cases$p
out <- t(sapply(seq_len(n), function(i) with(cases[i, ], {
  r <- encounter(E, pn, pt2, m_in, m_out, g, A = if (pn < 0) (m_in^2 + pt2)/(E - pn) else E + pn)
  c(match(r$type, c("miss", "reflect", "transmit")), r$E, r$pn, r$A)
})))
cases$type <- out[, 1]; cases$E_after <- out[, 2]; cases$pn_after <- out[, 3]; cases$A_after <- out[, 4]
write.csv(cases, "p1_cases.csv", row.names = FALSE)
cat(sprintf("\n3. %d random cases written for the 100-digit check (types: miss %d, reflect %d, transmit %d)\n",
            n, sum(out[, 1] == 1), sum(out[, 1] == 2), sum(out[, 1] == 3)))
# the 100-digit comparison of these cases (p1_reference.py), printed here so the letter can cite one script
ref <- system2("python3", c("p1_reference.py",
               "p1_cases.csv"), stdout = TRUE)
cat("\n4. AGAINST 100-DIGIT ARITHMETIC (p1_reference.py)\n"); cat(paste0("   ", ref), sep = "\n")
stopifnot(any(grepl("^PASS", ref)))
# how far below light speed a wall at gamma = 1e14 runs, quoted in the letter's method paragraph
cat(sprintf("\n5. AT gamma = 1e14 the wall's speed differs from light by 1 - beta = %.2g\n", 1/(1e28*(1 + sqrt(1 - 1e-28)))))
