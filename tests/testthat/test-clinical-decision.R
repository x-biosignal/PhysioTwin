# Clinical decision loop: personalise -> readiness gate -> score -> reconcile -> evaluate.

# one synthetic patient, personalised once and shared across the tests
.case <- local({
  truth <- limbTwin(strength = 6, damping = 0.6, stiffness = 3)
  imu <- imuMeasure(simulateTwin(truth, 5, 0.01),
                    imuSensor(gyro_noise = 0.03, accel_noise = 0.15), seed = 1)
  fit <- personalizeTwin(imu, limbTwin(strength = 4, damping = 0.6, stiffness = 3),
                         estimate = "strength", dt = 0.01)
  list(truth = truth, imu = imu, fit = fit, rdy = twinReadiness(fit, imu, dt = 0.01))
})

true_effect <- function(cand, th0)
  vapply(cand, function(c) {
    iv <- insilicoIntervention(.case$truth, scale = c$scale, set = c$set,
                               duration = 6, dt = 0.01, theta0 = th0)
    diff(range(iv$after$theta)) - diff(range(iv$before$theta))
  }, numeric(1))

test_that("interventions and the default library are well-formed", {
  skip_on_cran()
  iv <- intervention("x", scale = list(strength = 1.5), burden = 0.3, risk = 0.2)
  expect_s3_class(iv, "twin_intervention")
  lib <- defaultInterventions()
  expect_true(all(vapply(lib, inherits, logical(1), "twin_intervention")))
  expect_true("none" %in% names(lib))
  expect_output(print(iv), "Intervention")
})

test_that("the readiness gate passes an identifiable, well-fitting twin", {
  skip_on_cran()
  rdy <- .case$rdy
  expect_s3_class(rdy, "twin_readiness")
  expect_identical(rdy$verdict, "ready")
  expect_lt(rdy$reduced_chi2, 3)                              # fits within the noise
  expect_gt(rdy$confidence, 0.7)
  expect_lt(abs(.case$fit$estimates[["strength"]] - 6), 0.2) # recovered the truth
  expect_output(print(rdy), "READY")
})

test_that("the readiness gate rejects a twin that fits one channel but not the data", {
  skip_on_cran()
  # co-estimating two parameters from a poor start lands in a wrong local optimum that
  # reproduces the gyro but not the accelerometer -- the reduced chi-square exposes it
  bad <- personalizeTwin(.case$imu, limbTwin(strength = 4, stiffness = 2),
                         estimate = c("strength", "stiffness"), dt = 0.01)
  rdy <- twinReadiness(bad, .case$imu, dt = 0.01)
  expect_identical(rdy$verdict, "not_ready")
  expect_gt(rdy$reduced_chi2, 10)
  expect_lt(rdy$confidence, 0.1)
})

test_that("scoreInterventions recovers the mechanistic ordering with predictive intervals", {
  skip_on_cran()
  sc <- scoreInterventions(.case$fit$twin, defaultInterventions(), .case$rdy,
                           outcome = "rom", n_mc = 40, dt = 0.01, seed = 2)
  expect_s3_class(sc, "intervention_scores")
  r <- sc$ranking
  te <- true_effect(defaultInterventions(), .case$rdy$theta0)
  pe <- setNames(r$effect, r$intervention)[names(te)]
  expect_gt(cor(pe, te, method = "spearman"), 0.9)           # recovers the true ordering
  expect_equal(r$effect[r$intervention == "none"], 0, tolerance = 1e-6)
  expect_lt(r$effect[r$intervention == "orthosis"], 0)       # stiffening reduces ROM
  expect_true(all(r$ci_hi >= r$ci_lo))
  # a confident twin does not flag its (confident) recommendations, incl. confidently-ineffective ones
  expect_true(all(r$flag == ""))
  expect_output(print(sc), "Research prototype")             # honest-scope caveat is shown
})

test_that("evidence is reconciled and disagreement is flagged and penalised", {
  skip_on_cran()
  agree <- reconcileEvidence(0.30, 0.05, 0.27, 0.05)
  expect_identical(agree$agreement, "agree")
  expect_lt(agree$se, 0.05)                                  # combined is tighter than either input
  expect_gt(agree$estimate, 0.27); expect_lt(agree$estimate, 0.30)
  expect_identical(reconcileEvidence(0.30, 0.05, 0.05, 0.04)$agreement, "disagree")
  cand <- defaultInterventions()
  cand$strengthening$evidence <- list(effect = 0, se = 0.05) # contradicts the large mechanistic effect
  sc <- scoreInterventions(.case$fit$twin, cand, .case$rdy, n_mc = 30, dt = 0.01, seed = 3)
  row <- sc$ranking[sc$ranking$intervention == "strengthening", ]
  expect_identical(row$agree, "disagree")
  expect_identical(row$flag, "low-confidence")
})

test_that("the evaluation loop covers a realistic follow-up in its prediction interval", {
  skip_on_cran()
  te <- true_effect(list(defaultInterventions()$strengthening), .case$rdy$theta0)[1]
  ev <- evaluateIntervention(.case$fit$twin, defaultInterventions()$strengthening,
                             observed_delta = te + 0.05, .case$rdy,
                             outcome_sd = 0.1, n_mc = 60, dt = 0.01, seed = 6)
  expect_s3_class(ev, "intervention_evaluation")
  expect_true(ev$within_interval)                            # a plausible follow-up is covered
  expect_gt(ev$ci[2] - ev$ci[1], 0.1)                        # a usable (non-collapsed) interval
  ev2 <- evaluateIntervention(.case$fit$twin, defaultInterventions()$strengthening,
                              observed_delta = te + 5, .case$rdy, outcome_sd = 0.1,
                              n_mc = 60, dt = 0.01, seed = 6)
  expect_false(ev2$within_interval)                          # a wildly discordant outcome is caught
})

test_that("clinicalDecision runs the whole loop end to end", {
  skip_on_cran()
  dec <- clinicalDecision(.case$imu, limbTwin(strength = 4, damping = 0.6, stiffness = 3),
                          estimate = "strength", dt = 0.01, n_mc = 30, seed = 1)
  expect_s3_class(dec, "clinical_decision")
  expect_identical(dec$readiness$verdict, "ready")
  expect_identical(dec$scores$ranking$intervention[1], "strengthening")
  expect_output(print(dec), "Clinical decision loop")
})
