# Practical identifiability (profile likelihood) and Sobol sensitivity.

test_that("profile likelihood identifies a constrained parameter with a CI", {
  tw <- limbTwin(strength = 6, damping = 0.5, stiffness = 6)
  imu <- imuMeasure(simulateTwin(tw, 6, 0.03, theta0 = 0.2),
                    imuSensor(gyro_noise = 0.03, accel_noise = 0.12), seed = 1)
  pl <- profileLikelihood(imu, tw, "strength", nuisance = "damping",
                          over = c(0.5, 1.8), n_grid = 7, dt = 0.03)
  expect_s3_class(pl, "profile_likelihood")
  expect_true(pl$identifiable)                         # bounded profile -> identifiable
  expect_lt(abs(pl$best - 6) / 6, 0.15)               # best-fit near the true value
  expect_true(is.finite(pl$ci[1]) && is.finite(pl$ci[2]) && pl$ci[2] > pl$ci[1])
  expect_output(print(pl), "identifiable")
})

test_that("profile likelihood flags a structurally confounded parameter", {
  # g = 0, stiffness = 0 -> motion depends only on A/I and b/I, so strength,
  # inertia and damping are jointly unidentifiable
  tw <- limbTwin(strength = 6, damping = 0.4, stiffness = 0, g = 0)
  imu <- imuMeasure(simulateTwin(tw, 5, 0.03, theta0 = 0.2),
                    imuSensor(gyro_noise = 0.03, accel_noise = 0.12), seed = 2)
  pl <- profileLikelihood(imu, tw, "strength", nuisance = c("inertia", "damping"),
                          over = c(0.4, 2.5), n_grid = 9, dt = 0.03)
  expect_false(pl$identifiable)                        # flat profile -> not identifiable
  expect_lt(diff(range(pl$profile$cost)) / min(pl$profile$cost), 0.05)  # profile is flat
})

test_that("Sobol indices apportion output variance and flag a non-influential param", {
  so <- sobolIndices(list(strength = c(3, 9), damping = c(0.2, 1.0), stiffness = c(2, 6)),
                     n = 160, seed = 1)
  expect_s3_class(so, "sobol_indices")
  expect_true(all(so$Si >= 0 & so$Si <= 1))
  expect_true(all(so$STi >= 0 & so$STi <= 1))
  expect_true(all(so$STi >= so$Si - 0.15))             # total-effect >= first-order (estimator noise)
  # on range of motion, damping is influential, stiffness barely
  expect_gt(so$STi[so$param == "damping"], so$STi[so$param == "stiffness"])
  expect_output(print(so), "Sobol")
})

test_that("Sobol accepts a custom output metric", {
  so <- sobolIndices(list(strength = c(3, 9), damping = c(0.2, 1.0)),
                     output = function(tw) max(abs(simulateTwin(tw, 5, 0.03, theta0 = 0.2)$omega)),
                     n = 64, seed = 3)
  expect_equal(nrow(so), 2)
  expect_true(any(so$STi > 0.1))                       # something drives peak velocity
})
