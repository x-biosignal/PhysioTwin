# Layer 3: UKF data assimilation + twin personalisation.

test_that("UKF tracks a noisy random-walk signal", {
  skip_on_cran()
  set.seed(1)
  truth <- cumsum(rnorm(80, 0, 0.1))
  y <- matrix(truth + rnorm(80, 0, 0.2), 80, 1)
  fit <- unscentedKalmanFilter(function(x, k) x, function(x, k) x, y,
           Q = matrix(0.02), R = matrix(0.04), x0 = 0, P0 = matrix(1))
  expect_s3_class(fit, "ukf_result")
  # filtered estimate is closer to truth than the raw noisy measurement
  expect_lt(sqrt(mean((fit$state[, 1] - truth)^2)), sqrt(mean((y - truth)^2)))
})

test_that("optim personalisation recovers two parameters + gyro bias (robust)", {
  skip_on_cran()
  truth <- limbTwin(strength = 5.5, damping = 0.7, stiffness = 3.2)     # large-amplitude / nonlinear
  sim <- simulateTwin(truth, duration = 10, dt = 0.01, theta0 = 0.1)
  imu <- imuMeasure(sim, imuSensor(gyro_bias = 0.015, gyro_noise = 0.035,
                                   gyro_drift = 3e-4, accel_noise = 0.12), seed = 11)
  fit <- personalizeTwin(imu, limbTwin(strength = 4, damping = 0.3, stiffness = 3.2),
                         estimate = c("strength", "damping"), method = "optim", dt = 0.01,
                         gyro_var = 0.035^2, accel_var = 0.12^2)
  expect_s3_class(fit, "personalized_twin")
  expect_equal(unname(fit$estimates["strength"]), 5.5, tolerance = 0.5)
  expect_equal(unname(fit$estimates["damping"]),  0.7, tolerance = 0.2)
  expect_lt(abs(unname(fit$gyro_bias) - 0.015), 0.008)                  # nuisance bias recovered
  expect_output(print(fit), "Personalised movement twin")
})

test_that("UKF (recursive) recovers a parameter in the mild regime", {
  skip_on_cran()
  truth <- limbTwin(strength = 7, damping = 0.5, stiffness = 6)         # stiffer -> milder swing
  sim <- simulateTwin(truth, duration = 8, dt = 0.01, theta0 = 0.1)
  imu <- imuMeasure(sim, imuSensor(gyro_noise = 0.03, accel_noise = 0.12), seed = 1)
  fit <- personalizeTwin(imu, limbTwin(strength = 4, damping = 0.5, stiffness = 6),
                         estimate = "strength", method = "ukf", dt = 0.01,
                         gyro_var = 0.03^2, accel_var = 0.12^2)
  expect_equal(unname(fit$estimates["strength"]), 7, tolerance = 1.0)
  expect_false(is.null(fit$trajectory))                                # UKF returns a trajectory
})

test_that("the personalised twin reproduces the observed movement better than the prior", {
  skip_on_cran()
  truth <- limbTwin(strength = 8, damping = 0.4)
  sim <- simulateTwin(truth, duration = 8, dt = 0.01, theta0 = 0.1)
  imu <- imuMeasure(sim, imuSensor(gyro_noise = 0.03), seed = 5)
  prior <- limbTwin(strength = 4, damping = 0.4)
  fit <- personalizeTwin(imu, prior, estimate = "strength", dt = 0.01)
  rom_true <- diff(range(sim$theta))
  rom_prior <- diff(range(simulateTwin(prior, 8, 0.01, theta0 = 0.1)$theta))
  rom_fit  <- diff(range(simulateTwin(fit$twin, 8, 0.01, theta0 = 0.1)$theta))
  expect_lt(abs(rom_fit - rom_true), abs(rom_prior - rom_true))
})
