# Layer 4b: labelled synthetic training-data generation.

test_that("generateTrainingData returns features + labels in range", {
  d <- generateTrainingData(30, strength = c(3, 10), damping = c(0.1, 0.8),
                            duration = 4, dt = 0.02, seed = 1)
  expect_s3_class(d, "twin_dataset")
  expect_equal(nrow(d), 30)
  expect_true(all(c("gyro_rms", "ax_rms", "strength", "damping", "stiffness") %in% names(d)))
  expect_true(all(d$strength >= 3 & d$strength <= 10))
  expect_true(all(d$damping >= 0.1 & d$damping <= 0.8))
})

test_that("the synthetic data is learnable (features track the labels)", {
  d <- generateTrainingData(120, duration = 5, dt = 0.02, randomize_sensor = TRUE, seed = 2)
  # gyro amplitude should rise with muscle strength despite sensor randomisation
  expect_gt(cor(d$gyro_rms, d$strength), 0.4)
  expect_output(print(d), "Twin dataset")
})

test_that("raw-stream mode returns labels + IMU streams", {
  d <- generateTrainingData(5, duration = 3, dt = 0.02, features = FALSE, seed = 3)
  expect_equal(nrow(d$labels), 5)
  expect_length(d$streams, 5)
  expect_true(all(c("gyro", "ax", "ay") %in% names(d$streams[[1]])))
})
