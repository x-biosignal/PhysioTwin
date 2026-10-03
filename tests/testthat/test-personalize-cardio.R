# Multi-system personalisation: fit the closed-loop CV-respiratory twin to HRV/BPV.

test_that("cardioSummary extracts the HRV/BPV summaries", {
  skip_on_cran()
  s <- cardioSummary(cardioRespiratory(duration = 120, seed = 1))
  expect_named(s, c("brs", "bp_lf", "hr_hf"))
  expect_true(all(is.finite(s)))
})

test_that("personalizeCardio fits the closed loop to reproduce a subject's HRV/BPV", {
  skip_on_cran()
  truep <- list(g_vagal = 0.013, g_resist = 0.008, rsa_amp = 0.045)
  # the subject's summaries (averaged over the seeds the fit will use, so the target is stable)
  S <- vapply(2:4, function(s) cardioSummary(do.call(cardioRespiratory,
              c(truep, list(duration = 300, seed = s)))), numeric(3))
  target <- rowMeans(S); names(target) <- c("brs", "bp_lf", "hr_hf")
  fit <- personalizeCardio(target, n_init = 10, n_iter = 14, n_rep = 3, seed = 1)
  expect_s3_class(fit, "cardio_personalization")
  # the fitted twin reproduces the subject's baroreflex sensitivity and LF/HF variability
  expect_lt(max(abs(fit$rel_error)), 0.15)
  expect_true(all(names(fit$estimates) %in% c("g_vagal", "g_resist", "rsa_amp")))
  expect_output(print(fit), "personalisation")
})
