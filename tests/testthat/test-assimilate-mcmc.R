# Bayesian posterior engines: adaptive Metropolis-Hastings and sequential ABC.

test_that("metropolis samples a correlated bivariate-normal target", {
  mu <- c(1, -2); Sig <- matrix(c(1, 0.6, 0.6, 1.5), 2, 2); Si <- solve(Sig)
  lp <- function(x) { d <- x - mu; as.numeric(-0.5 * t(d) %*% Si %*% d) }
  m <- metropolis(lp, init = c(a = 0, b = 0), n_iter = 8000, seed = 1)
  expect_s3_class(m, "mcmc_result")
  expect_lt(max(abs(m$mean - mu)), 0.15)                # recovers the mean ...
  expect_lt(abs(var(m$samples[, 1]) - 1), 0.3)          # ... the variances ...
  expect_lt(abs(var(m$samples[, 2]) - 1.5), 0.4)
  expect_lt(abs(cov(m$samples[, 1], m$samples[, 2]) - 0.6), 0.3)   # ... and the covariance
  expect_gt(m$acceptance, 0.1); expect_lt(m$acceptance, 0.6)
  expect_output(print(m), "MCMC")
})

test_that("metropolis recovers the analytic normal-mean posterior", {
  set.seed(2); y <- rnorm(30, 3, 1)
  pv <- 1 / (30 + 1 / 100); pm <- pv * sum(y)           # N(0,100) prior, unit variance
  lp <- function(mm) sum(dnorm(y, mm, 1, log = TRUE)) + dnorm(mm, 0, 10, log = TRUE)
  m <- metropolis(lp, init = c(mu = 0), n_iter = 6000, seed = 3)
  expect_lt(abs(unname(m$mean) - pm), 0.1)
  expect_lt(abs(unname(m$sd) - sqrt(pv)), 0.04)
})

test_that("SMC-ABC recovers a Gaussian with a shrinking tolerance", {
  set.seed(4); dat <- rnorm(300, 2, 1.5); obs <- c(mean(dat), sd(dat))
  sim <- function(th) { x <- rnorm(300, th[1], th[2]); c(mean(x), sd(x)) }
  pri <- function() c(m = runif(1, -3, 7), s = runif(1, 0.2, 4))
  sm <- abcSMC(sim, obs, pri, n_particles = 150, n_populations = 6, quantile = 0.4, seed = 1)
  expect_s3_class(sm, "abc_smc")
  expect_lt(abs(sm$mean[["m"]] - obs[1]), 0.3)          # recovers the observed summaries
  expect_lt(abs(sm$mean[["s"]] - obs[2]), 0.3)
  expect_lt(sm$tolerance[6], sm$tolerance[1])           # the tolerance schedule shrinks
  expect_lt(sm$sd[["m"]], 0.5 * (10 / sqrt(12)))        # posterior concentrates far below the prior
  expect_output(print(sm), "SMC-ABC")
})
