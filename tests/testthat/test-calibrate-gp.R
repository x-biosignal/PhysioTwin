# GP-accelerated calibration (Bayesian optimisation).

test_that("gpCalibrate finds a 2-D minimum in few evaluations and beats random search", {
  skip_on_cran()
  loss <- function(p) (p[["x"]] - 0.3)^2 + (p[["y"]] - 0.7)^2
  cal <- gpCalibrate(loss, list(x = c(0, 1), y = c(0, 1)), n_init = 8, n_iter = 18, seed = 1)
  expect_s3_class(cal, "gp_calibration")
  expect_lt(sqrt((cal$par[["x"]] - 0.3)^2 + (cal$par[["y"]] - 0.7)^2), 0.1)   # near the true optimum
  set.seed(9); rnd <- min(apply(matrix(runif(26 * 2), 26, 2), 1,
                                function(r) (r[1] - 0.3)^2 + (r[2] - 0.7)^2))
  expect_lt(cal$value, rnd)                              # better than equal-budget random search
  expect_lte(cal$trace[18], cal$trace[1])               # best-so-far is non-increasing
  expect_output(print(cal), "GP calibration")
})

test_that("gpCalibrate calibrates a twin parameter to a target outcome", {
  skip_on_cran()
  target <- diff(range(simulateTwin(limbTwin(strength = 6.3), duration = 6, dt = 0.01)$theta))
  loss <- function(p) {
    r <- diff(range(simulateTwin(limbTwin(strength = p[["strength"]]), duration = 6, dt = 0.01)$theta))
    (r - target)^2
  }
  cal <- gpCalibrate(loss, list(strength = c(3, 9)), n_init = 6, n_iter = 12, seed = 2)
  expect_lt(abs(cal$par[["strength"]] - 6.3), 0.3)      # recovers the true parameter
})

test_that("the GP emulator tolerates near-duplicate design points", {
  skip_on_cran()
  set.seed(3); X <- rbind(matrix(runif(20), 10, 2), c(0.5, 0.5), c(0.5 + 1e-9, 0.5))  # a near-duplicate
  y <- apply(X, 1, function(r) sum((r - 0.4)^2))
  gp <- gpEmulator(X, y)                                 # must not fail on the singular kernel
  expect_s3_class(gp, "gpEmulator")
  expect_true(all(is.finite(predict(gp, X)$mean)))
})
