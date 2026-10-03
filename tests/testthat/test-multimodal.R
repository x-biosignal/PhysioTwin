# Multi-modal coupling: one latent driver -> movement + cardiorespiratory.

test_that("fitness and effort map coherently to movement and physiology", {
  lo <- multimodalTwin(fitness = 0.2, effort = 0.5)
  hi <- multimodalTwin(fitness = 0.8, effort = 0.5)
  expect_s3_class(hi, "multimodal_twin")
  expect_gt(hi$movement$strength, lo$movement$strength)     # fitter -> stronger
  expect_lt(hi$movement$damping, lo$movement$damping)       # fitter -> less damped
  expect_lt(hi$resting_hr, lo$resting_hr)                   # fitter -> lower resting HR
  # effort raises heart rate and breathing
  expect_gt(multimodalTwin(0.5, effort = 0.9)$resting_hr, multimodalTwin(0.5, effort = 0.1)$resting_hr)
  expect_gt(multimodalTwin(0.5, effort = 0.9)$resp_rate, multimodalTwin(0.5, effort = 0.1)$resp_rate)
  expect_output(print(hi), "Multi-modal twin")
})

test_that("simulateMultimodal returns a coherent cross-modal summary", {
  s <- simulateMultimodal(multimodalTwin(fitness = 0.6), duration_hr = 60, seed = 1)
  expect_named(s, c("rom", "peak_velocity", "mean_hr", "resp_rate"))
  expect_true(all(is.finite(s)))
  expect_equal(unname(s["mean_hr"]), multimodalTwin(fitness = 0.6)$resting_hr, tolerance = 3)
})

test_that("a fitness intervention moves movement and heart rate in the right directions", {
  mt <- multimodalTwin(fitness = 0.4)
  iv <- insilicoInterventionMultimodal(mt, fitness_gain = 0.3, duration_hr = 60, seed = 1)
  expect_s3_class(iv, "mm_intervention")
  expect_gt(iv$change["rom"], 0)                            # rehab raises range of motion
  expect_lt(iv$change["mean_hr"], 0)                        # ... and lowers resting heart rate
  expect_output(print(iv), "Multi-modal intervention")
})

test_that("the cross-modal effect grows with the intervention size", {
  mt <- multimodalTwin(fitness = 0.3)
  small <- insilicoInterventionMultimodal(mt, fitness_gain = 0.2, duration_hr = 60, seed = 1)
  large <- insilicoInterventionMultimodal(mt, fitness_gain = 0.5, duration_hr = 60, seed = 1)
  expect_gt(large$change["rom"], small$change["rom"])
  expect_lt(large$change["mean_hr"], small$change["mean_hr"])
})
