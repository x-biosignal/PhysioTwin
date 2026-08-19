# personalizeTwin multistart: robust multi-parameter convergence from a poor prior.

test_that("multistart recovers a two-parameter fit that a single start misses", {
  skip_on_cran()
  truth <- limbTwin(strength = 6, damping = 0.6, stiffness = 3)
  imu <- imuMeasure(simulateTwin(truth, 5, 0.02),
                    imuSensor(gyro_noise = 0.03, accel_noise = 0.1), seed = 1)
  poor <- limbTwin(strength = 12, damping = 0.6, stiffness = 1)    # a deliberately poor guess
  single <- personalizeTwin(imu, poor, estimate = c("strength", "stiffness"),
                            dt = 0.02, n_starts = 1)
  multi  <- personalizeTwin(imu, poor, estimate = c("strength", "stiffness"), dt = 0.02)
  # the grid pre-scan lands the fit in the true basin; a lone start gets stuck worse
  expect_lt(abs(multi$estimates[["strength"]] - 6), 0.5)
  expect_lt(abs(multi$estimates[["stiffness"]] - 3), 0.6)
  expect_lt(multi$objective, single$objective)
})

test_that("multistart keeps single-parameter recovery and is deterministic", {
  skip_on_cran()
  truth <- limbTwin(strength = 6, damping = 0.6, stiffness = 3)
  imu <- imuMeasure(simulateTwin(truth, 5, 0.02),
                    imuSensor(gyro_noise = 0.03, accel_noise = 0.1), seed = 1)
  prior <- limbTwin(strength = 4, damping = 0.6, stiffness = 3)
  fit  <- personalizeTwin(imu, prior, estimate = "strength", dt = 0.02)
  fit2 <- personalizeTwin(imu, prior, estimate = "strength", dt = 0.02)
  expect_lt(abs(fit$estimates[["strength"]] - 6), 0.3)
  expect_identical(fit$estimates, fit2$estimates)             # deterministic grid -> reproducible
})
