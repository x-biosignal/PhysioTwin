# Whole-record joint fitting of the closed-loop cardiovascular state-space model.

test_that("cardioStateSpace builds a valid linear-Gaussian state-space model", {
  skip_on_cran()
  resp <- sin(2 * pi * 0.25 * cumsum(rep(0.85, 100)))
  m <- cardioStateSpace(g = 8, rho = 0.9, A_rsa = 20, resp = resp)
  expect_true(is.function(m$f) && is.function(m$h))
  expect_equal(dim(m$Q), c(2, 2)); expect_length(m$x0, 2)
  expect_equal(m$h(c(1, 0), 1), 850 + 8 * 1 + 20 * resp[1])   # rr0 + g*state + A_rsa*resp
})

test_that("whole-record fitting jointly recovers the closed-loop parameters", {
  skip_on_cran()
  set.seed(2); T <- 200; rr0 <- 850; g <- 8; rho <- 0.9; A_rsa <- 20; om <- 0.53
  A <- matrix(c(2 * rho * cos(om), 1, -rho^2, 0), 2, 2)
  x <- c(0, 0); rrm <- numeric(T)                          # Mayer-driven RR ...
  for (k in seq_len(T)) { x <- A %*% x + c(rnorm(1, 0, 1), 0); rrm[k] <- rr0 + g * x[1] + rnorm(1, 0, 3) }
  resp <- sin(2 * pi * 0.25 * (0:(T - 1)) * mean(rrm) / 1000)
  rr <- rrm + A_rsa * resp                                 # ... plus the respiratory (RSA) term
  fit <- personalizeCardioWaveform(rr, resp_freq = 0.25, estimate = c("g", "rho", "A_rsa"),
                                   fixed = list(omega_m = om, rr0 = NULL, q_mayer = 1, r_obs = 9),
                                   n_iter = 2500, seed = 1)
  expect_s3_class(fit, "cardio_waveform_fit")
  q <- fit$posterior$quantiles
  # the whole record identifies all three parameters -- where the summary fit had a ridge
  expect_true(8 > q[1, 1] && 8 < q[3, 1])                  # Mayer gain g
  expect_true(0.9 > q[1, 2] && 0.9 < q[3, 2])              # Mayer damping rho
  expect_true(20 > q[1, 3] && 20 < q[3, 3])                # RSA amplitude
  expect_gt(fit$posterior$acceptance, 0.1)
  expect_output(print(fit), "whole-record")
})

test_that("the broadband (ar1) low-frequency model recovers a peakless LF record", {
  skip_on_cran()
  # RR with BROADBAND low-frequency variability (a first-order AR, no Mayer peak) + RSA
  set.seed(3); T <- 300; rr0 <- 850; phi <- 0.9; g <- 6; A_rsa <- 15
  x <- 0; rrm <- numeric(T)
  for (k in seq_len(T)) { x <- phi * x + rnorm(1, 0, 1); rrm[k] <- rr0 + g * x + rnorm(1, 0, 3) }
  resp <- sin(2 * pi * 0.25 * (0:(T - 1)) * mean(rrm) / 1000); rr <- rrm + A_rsa * resp
  m <- cardioStateSpace(g = g, A_rsa = A_rsa, resp = resp, lf_model = "ar1", phi = phi)
  expect_identical(m$lf_model, "ar1"); expect_length(m$x0, 1)          # first-order (1-D) state
  fit <- personalizeCardioWaveform(rr, resp_freq = 0.25, lf_model = "ar1",
                                   fixed = list(rr0 = NULL, q_mayer = 1, r_obs = 9),
                                   n_iter = 3000, seed = 1)
  expect_identical(fit$lf_model, "ar1")
  expect_identical(fit$estimate, c("g", "phi", "A_rsa"))              # no arbitrary Mayer frequency
  expect_lt(abs(fit$estimates[["g"]] - g), 1.5)
  expect_lt(abs(fit$estimates[["phi"]] - phi), 0.06)
  expect_lt(abs(fit$estimates[["A_rsa"]] - A_rsa), 3)
})
