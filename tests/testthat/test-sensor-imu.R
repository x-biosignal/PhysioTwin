# Layer 2: the forward IMU sensor model.

test_that("noise-free gyro equals the true angular velocity", {
  sim <- simulateTwin(limbTwin(), duration = 3, dt = 0.01)
  imu <- imuMeasure(sim, imuSensor(), noise = FALSE)
  expect_equal(imu$gyro, sim$omega, tolerance = 1e-12)
})

test_that("a static sensor reads gravity (accel magnitude = g)", {
  f <- PhysioTwin:::.imu_forward
  a <- f(theta = 0.3, omega = 0, alpha = 0, r = 0.25, g = 9.81)  # static: only gravity
  expect_equal(sqrt(a["ax"]^2 + a["ay"]^2), 9.81, tolerance = 1e-9, ignore_attr = TRUE)
  # horizontal, no motion -> zero specific force
  a0 <- f(theta = 0.2, omega = 0, alpha = 0, r = 0.25, g = 0)
  expect_equal(unname(a0[c("ax", "ay")]), c(0, 0), tolerance = 1e-12)
})

test_that("sensor imperfections add bias and noise", {
  sim <- simulateTwin(limbTwin(), duration = 4, dt = 0.01)
  clean <- imuMeasure(sim, noise = FALSE)
  noisy <- imuMeasure(sim, imuSensor(gyro_bias = 0.1, gyro_noise = 0.05), seed = 1)
  expect_gt(mean(noisy$gyro - clean$gyro), 0.05)       # positive bias shifts the gyro
  expect_gt(sd(noisy$gyro - clean$gyro), 0.01)         # noise present
  expect_true(all(c("time", "gyro", "ax", "ay") %in% names(noisy)))
})

test_that("random-walk drift accumulates over time", {
  sim <- simulateTwin(limbTwin(), duration = 10, dt = 0.01)
  d <- imuMeasure(sim, imuSensor(gyro_noise = 0, gyro_drift = 2e-3), seed = 2)
  clean <- imuMeasure(sim, noise = FALSE)
  dev <- d$gyro - clean$gyro
  expect_gt(sd(dev[900:1000]), sd(dev[1:100]))         # drift grows
  expect_output(print(imuSensor()), "IMU sensor")
})
