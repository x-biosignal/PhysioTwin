# Gaussian-process emulator (surrogate model).

test_that("the GP interpolates its training runs and predicts a held-out function", {
  set.seed(3); f <- function(x) sin(6 * x) + 0.3 * x
  Xtr <- matrix(sort(runif(16)), 16, 1); ytr <- f(Xtr[, 1])
  gp <- gpEmulator(Xtr, ytr)
  expect_s3_class(gp, "gpEmulator")
  tr <- predict(gp, Xtr)
  expect_named(tr, c("mean", "sd"))
  expect_lt(max(abs(tr$mean - ytr)), 1e-3)                    # near-exact interpolation
  Xte <- seq(0.03, 0.97, length.out = 120)
  expect_lt(sqrt(mean((predict(gp, Xte)$mean - f(Xte))^2)), 0.02)   # accurate emulation
  expect_output(print(gp), "GP emulator")
})

test_that("predictive uncertainty grows away from the training data", {
  set.seed(3); f <- function(x) sin(6 * x) + 0.3 * x
  Xtr <- matrix(sort(runif(16)), 16, 1); gp <- gpEmulator(Xtr, f(Xtr[, 1]))
  sd_in  <- predict(gp, 0.5)$sd                               # inside the training range
  sd_out <- predict(gp, 2.0)$sd                               # far extrapolation
  expect_gt(sd_out, 20 * sd_in)
  expect_lt(abs(sd_out - gp$sigma_f), 0.3 * gp$sigma_f)       # approaches the prior signal SD
})

test_that("a fixed nugget gives exact interpolation and vector input works", {
  set.seed(4); x <- runif(12, 0, 6); y <- cos(x)
  gp <- gpEmulator(x, y, optimize = FALSE, sigma_n = 1e-6)    # vector input, tiny noise
  p <- predict(gp, x)
  expect_length(p$mean, 12)
  expect_lt(max(abs(p$mean - y)), 1e-3)
})
