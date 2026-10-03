# The validity-assessment harness (verification & validation).

test_that("L0 verification: energy conserved and RK4 is ~4th order", {
  skip_on_cran()
  v <- validateTwin(limbTwin(strength = 6, damping = 0.5, stiffness = 6),
                    checks = "verification")
  expect_lt(v$verification$energy_drift, 1e-3)         # dissipation-free -> conserved
  expect_gt(v$verification$rk4_order, 3.5)             # ~4 for RK4
  expect_lt(v$verification$rk4_order, 4.5)
})

test_that("L1 identifiability: standard parameters are structurally identifiable", {
  skip_on_cran()
  v <- validateTwin(limbTwin(strength = 6, damping = 0.5, stiffness = 6),
                    estimate = c("strength", "damping"), checks = "identifiability")
  expect_equal(v$identifiability$rank, 2)              # full rank -> identifiable
  expect_true(is.finite(v$identifiability$condition))
  expect_length(v$identifiability$cramer_rao_sd, 2)
})

test_that("L2 recovery: optim assimilation is unbiased with high coverage", {
  skip_on_cran()
  v <- validateTwin(limbTwin(strength = 6, damping = 0.5, stiffness = 6),
                    estimate = c("strength", "damping"), method = "optim",
                    n_mc = 10, checks = "recovery")
  expect_lt(abs(v$recovery$bias[v$recovery$param == "strength"]), 0.5)
  expect_gt(v$recovery$coverage_15pct[v$recovery$param == "strength"], 0.7)
})

test_that("L3 sensitivity: influential parameters get positive elementary effects", {
  skip_on_cran()
  v <- validateTwin(limbTwin(strength = 6, damping = 0.5, stiffness = 6),
                    estimate = c("strength", "damping"), n_traj = 6, checks = "sensitivity")
  expect_true(all(v$sensitivity >= 0))
  expect_true(any(v$sensitivity > 0.1))               # at least one clearly matters
})

test_that("L4 prediction: the fitted twin predicts an unseen condition + face validity", {
  skip_on_cran()
  v <- validateTwin(limbTwin(strength = 6, damping = 0.5, stiffness = 6),
                    estimate = c("strength", "damping"), checks = "prediction")
  expect_lt(v$prediction$cross_condition_rel_rmse, 0.2)   # generalises to a new frequency
  expect_true(all(v$prediction$face_validity))            # known monotone behaviours hold
})

test_that("L5 domain: amplitude and regime are reported; full report prints", {
  skip_on_cran()
  v <- validateTwin(limbTwin(strength = 6, damping = 0.5, stiffness = 6),
                    estimate = c("strength", "damping"), n_mc = 6, n_traj = 4)
  expect_s3_class(v, "twin_validity")
  expect_true(v$domain$amplitude_rad > 0)
  expect_true(v$domain$regime %in% c("moderate", "large-amplitude (strongly nonlinear)"))
  expect_output(print(v), "Twin validity report")
})
