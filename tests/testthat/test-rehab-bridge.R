# PhysioRehab bridge: hold the twin's prediction to a real clinical SCED/MCID/ICF evaluation.

test_that("rehabEvaluate runs the clinical SCED/MCID/ICF evaluation and reconciles the twin", {
  skip_if_not_installed("PhysioRehab")
  set.seed(1)
  base <- rnorm(5, 0.80, 0.03); post <- rnorm(6, 1.00, 0.03)     # a clear improvement
  ev <- rehabEvaluate(base, post, predicted_delta = 0.18, measure = "gait_speed", mcid = 0.16)
  expect_s3_class(ev, "rehab_evaluation")
  expect_gt(ev$sced$nap, 0.9)                        # strong single-case effect
  expect_true(isTRUE(ev$sced$exceeds_mcid))          # clinically meaningful change
  expect_identical(ev$icf, "d450")                   # gait speed -> ICF activity category
  expect_true(ev$reconciliation$sign_agree)
  expect_true(ev$reconciliation$mcid_agree)          # twin and reality agree it exceeds MCID
  expect_output(print(ev), "SCED")
})

test_that("rehabEvaluate flags a twin over-prediction against the real series", {
  skip_if_not_installed("PhysioRehab")
  set.seed(2)
  base <- rnorm(5, 0.80, 0.03); post <- rnorm(6, 0.86, 0.03)     # a real but sub-MCID change
  ev <- rehabEvaluate(base, post, predicted_delta = 0.18, measure = "gait_speed", mcid = 0.16)
  expect_false(isTRUE(ev$sced$exceeds_mcid))         # reality: below the MCID
  expect_true(ev$reconciliation$sign_agree)          # the twin had the direction right ...
  expect_false(ev$reconciliation$mcid_agree)         # ... but over-predicted the magnitude
  expect_lt(ev$reconciliation$observed, ev$reconciliation$predicted)
})

test_that("rehabEvaluate returns the clinical result alone when no twin prediction is given", {
  skip_if_not_installed("PhysioRehab")
  ev <- rehabEvaluate(rnorm(4, 1, 0.05), rnorm(5, 1.3, 0.05), measure = "gait_speed")
  expect_null(ev$reconciliation)
  expect_true(is.list(ev$sced))
  expect_true(is.numeric(ev$sced$delta))
})
