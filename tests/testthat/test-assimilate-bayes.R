# Approximate Bayesian Computation: likelihood-free parameter calibration.

test_that("ABC recovers the parameters of a Gaussian and covers the truth", {
  set.seed(1); data <- rnorm(60, mean = 3, sd = 1); obs <- c(mean(data), sd(data))
  sim <- function(th) { x <- rnorm(60, th["mu"], th["sig"]); c(mean(x), sd(x)) }
  pri <- function() c(mu = runif(1, -3, 9), sig = runif(1, 0.2, 3))
  post <- abcCalibration(sim, obs, pri, n = 6000, accept = 0.03, seed = 2)
  expect_s3_class(post, "abc_posterior")
  expect_equal(nrow(post$samples), 180)                       # 3% of 6000
  expect_lt(abs(post$mean["mu"] - 3), 0.3)                    # posterior mean near truth
  expect_lt(abs(post$mean["sig"] - 1), 0.3)
  # 95% credible interval covers the truth
  expect_lt(post$quantiles[1, 1], 3); expect_gt(post$quantiles[3, 1], 3)
  expect_output(print(post), "ABC posterior")
})

test_that("a tighter tolerance concentrates the posterior", {
  set.seed(1); obs <- c(mean(rnorm(60, 3, 1)), 1)
  sim <- function(th) { x <- rnorm(60, th["mu"], th["sig"]); c(mean(x), sd(x)) }
  pri <- function() c(mu = runif(1, -3, 9), sig = runif(1, 0.2, 3))
  loose <- abcCalibration(sim, obs, pri, n = 6000, accept = 0.05, seed = 2)
  tight <- abcCalibration(sim, obs, pri, n = 6000, accept = 0.01, seed = 2)
  expect_lt(tight$sd["mu"], loose$sd["mu"])
})

test_that("the regression adjustment improves on plain rejection", {
  set.seed(1); data <- rnorm(60, 3, 1); obs <- c(mean(data), sd(data))
  sim <- function(th) { x <- rnorm(60, th["mu"], th["sig"]); c(mean(x), sd(x)) }
  pri <- function() c(mu = runif(1, -3, 9), sig = runif(1, 0.2, 3))
  reg   <- abcCalibration(sim, obs, pri, n = 6000, accept = 0.05, regression = TRUE,  seed = 2)
  plain <- abcCalibration(sim, obs, pri, n = 6000, accept = 0.05, regression = FALSE, seed = 2)
  expect_lt(abs(reg$mean["mu"] - 3), abs(mean(plain$raw[, "mu"]) - 3) + 1e-9)
})

test_that("abcCalibration is reproducible under a fixed seed", {
  sim <- function(th) { x <- rnorm(50, th["mu"], 1); c(mean(x), sd(x)) }
  pri <- function() c(mu = runif(1, -3, 9))
  expect_identical(abcCalibration(sim, c(2, 1), pri, n = 1000, seed = 5)$samples,
                   abcCalibration(sim, c(2, 1), pri, n = 1000, seed = 5)$samples)
})
