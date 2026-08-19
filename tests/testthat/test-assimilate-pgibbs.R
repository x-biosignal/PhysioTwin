# Particle Gibbs with ancestor sampling: joint state-and-parameter sampling.

kf_loglik_ar1 <- function(a, y, q, r) {
  x <- 0; P <- 1; ll <- 0
  for (k in seq_along(y)) {
    xp <- a * x; Pp <- a^2 * P + q; S <- Pp + r; v <- y[k] - xp
    ll <- ll - 0.5 * (v^2 / S + log(2 * pi * S)); K <- Pp / S; x <- xp + K * v; P <- (1 - K) * Pp
  }
  ll
}

test_that("particle Gibbs recovers the posterior and the latent trajectory", {
  skip_on_cran()
  set.seed(1); a <- 0.8; xt <- numeric(45); for (k in 2:45) xt[k] <- a * xt[k - 1] + rnorm(1, 0, 0.3)
  y <- xt + rnorm(45, 0, 0.3); q <- 0.09; r <- 0.09
  # exact parameter posterior via the Kalman likelihood
  ag <- seq(-0.98, 0.98, 0.01)
  lp <- vapply(ag, function(aa) kf_loglik_ar1(aa, y, q, r) + dunif(aa, -1, 1, log = TRUE), numeric(1))
  w <- exp(lp - max(lp)); w <- w / sum(w); em <- sum(ag * w); es <- sqrt(sum(w * (ag - em)^2))
  bm <- function(th) list(f = function(x, k) th[1] * x, h = function(x, k) x,
                          Q = matrix(q), R = matrix(r), x0 = 0, P0 = matrix(1))
  fit <- particleGibbs(bm, matrix(y, 45, 1), function(th) dunif(th[1], -1, 1, log = TRUE),
                       init = c(a = 0.3), n_particles = 50, n_iter = 600, seed = 2)
  expect_s3_class(fit, "pgibbs_result")
  expect_lt(abs(fit$mean[["a"]] - em), 0.1)              # matches the exact-likelihood posterior
  expect_lt(abs(fit$sd[["a"]] - es), 0.06)
  # the conditional-SMC sweep also returns a state trajectory that tracks the truth
  expect_equal(dim(fit$trajectory), c(45, 1))
  expect_gt(cor(fit$trajectory[, 1], xt), 0.75)
  expect_gt(fit$acceptance, 0.1); expect_lt(fit$acceptance, 0.7)
  expect_output(print(fit), "Particle Gibbs")
})
