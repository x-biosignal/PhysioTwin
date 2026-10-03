# Particle MCMC (PMMH): whole-waveform Bayesian parameter estimation.

# exact Kalman marginal log-likelihood of an AR(1)+Gaussian state-space model
kf_loglik <- function(a, y, q, r) {
  x <- 0; P <- 1; ll <- 0
  for (k in seq_along(y)) {
    xp <- a * x; Pp <- a^2 * P + q; S <- Pp + r; v <- y[k] - xp
    ll <- ll - 0.5 * (v^2 / S + log(2 * pi * S)); K <- Pp / S; x <- xp + K * v; P <- (1 - K) * Pp
  }
  ll
}

test_that("particle MCMC matches the exact-likelihood posterior on a linear-Gaussian model", {
  skip_on_cran()
  set.seed(1); a <- 0.7; x <- numeric(50); for (k in 2:50) x[k] <- a * x[k - 1] + rnorm(1, 0, 0.3)
  y <- x + rnorm(50, 0, 0.3); q <- 0.09; r <- 0.09
  # exact posterior of the AR coefficient via the Kalman likelihood on a grid
  ag <- seq(-0.98, 0.98, 0.01)
  lp <- vapply(ag, function(aa) kf_loglik(aa, y, q, r) + dunif(aa, -1, 1, log = TRUE), numeric(1))
  w <- exp(lp - max(lp)); w <- w / sum(w)
  em <- sum(ag * w); es <- sqrt(sum(w * (ag - em)^2))
  bm <- function(th) list(f = function(x, k) th[1] * x, h = function(x, k) x,
                          Q = matrix(q), R = matrix(r), x0 = 0, P0 = matrix(1))
  fit <- particleMCMC(bm, matrix(y, 50, 1), function(th) dunif(th[1], -1, 1, log = TRUE),
                      init = c(a = 0.3), n_particles = 100, n_iter = 900, seed = 2)
  expect_s3_class(fit, "pmcmc_result")
  expect_lt(abs(fit$mean[["a"]] - em), 0.08)         # particle-filter-likelihood posterior ...
  expect_lt(abs(fit$sd[["a"]] - es), 0.05)           # ... matches the exact-likelihood one
  expect_gt(fit$acceptance, 0.1); expect_lt(fit$acceptance, 0.7)
  expect_output(print(fit), "Particle MCMC")
})

test_that("particle MCMC recovers a baroreflex gain from the whole RR waveform", {
  skip_on_cran()
  set.seed(5); phi <- 0.6; g <- 0.8; s <- numeric(70); for (k in 2:70) s[k] <- phi * s[k - 1] + rnorm(1, 0, 1)
  rr <- g * s + rnorm(70, 0, 0.3)                    # observed RR deviation = gain x latent pressure
  bm <- function(th) list(f = function(x, k) phi * x, h = function(x, k) th[1] * x,
                          Q = matrix(1), R = matrix(0.09), x0 = 0, P0 = matrix(1 / (1 - phi^2)))
  fit <- particleMCMC(bm, matrix(rr, 70, 1), function(th) dunif(th[1], 0, 3, log = TRUE),
                      init = c(gain = 0.4), n_particles = 100, n_iter = 700, seed = 3)
  expect_lt(abs(fit$mean[["gain"]] - 0.8), 0.25)     # recovered from the whole series ...
  expect_gt(0.8, fit$quantiles[1, 1]); expect_lt(0.8, fit$quantiles[3, 1])   # ... truth in the 95% CI
})
