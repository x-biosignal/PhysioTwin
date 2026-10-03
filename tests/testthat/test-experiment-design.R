# Bayesian optimal experimental design: choose the most informative condition.

test_that("optimalDesign ranks driving frequencies by Fisher information", {
  skip_on_cran()
  od <- optimalDesign(limbTwin(strength = 6), estimate = "strength",
                      designs = seq(0.4, 2.4, by = 0.2), design_var = "freq",
                      sensor = imuSensor(gyro_noise = 0.05, accel_noise = 0.2),
                      duration = 5, dt = 0.02)
  expect_s3_class(od, "optimal_design")
  expect_equal(nrow(od$table), 11)
  expect_true(all(is.finite(od$table$logdet_fim)))
  expect_gt(max(od$table$logdet_fim) - min(od$table$logdet_fim), 3)   # non-degenerate information curve
  expect_equal(od$best, od$table$design[which.max(od$table$score)])
  expect_lt(od$best, 1.2)                                             # optimum near the twin's resonance (~0.65 Hz)
  expect_equal(od$table$crlb_strength[od$table$design == od$best],
               min(od$table$crlb_strength))                          # smallest Cramer-Rao SD
  expect_output(print(od), "Optimal experimental design")
})

test_that("the optimal design's Cramer-Rao bound predicts real recovery variance", {
  skip_on_cran()
  truth <- limbTwin(strength = 6); sens <- imuSensor(gyro_noise = 0.05, accel_noise = 0.2)
  od <- optimalDesign(truth, estimate = "strength", designs = c(0.6, 2.4),
                      design_var = "freq", sensor = sens, duration = 4, dt = 0.02)
  expect_lt(od$table$crlb_strength[od$table$design == 0.6],
            od$table$crlb_strength[od$table$design == 2.4])           # CRLB: 0.6 Hz far more informative
  # ... and personalising there really recovers strength with lower empirical variance
  rec_sd <- function(fr) {
    tw <- truth; tw$freq <- fr
    sd(vapply(1:8, function(i) {
      imu <- imuMeasure(simulateTwin(tw, 4, 0.02, theta0 = 0.1), sens, seed = 100 + i)
      personalizeTwin(imu, limbTwin(strength = 4, freq = fr), estimate = "strength",
                      dt = 0.02, n_starts = 1, gyro_var = 0.05^2, accel_var = 0.2^2)$estimates[["strength"]]
    }, numeric(1)))
  }
  expect_lt(rec_sd(0.6), rec_sd(2.4))
})
