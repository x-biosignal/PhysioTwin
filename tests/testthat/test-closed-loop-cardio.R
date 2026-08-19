# Closed-loop cardiovascular-respiratory twin: emergent Mayer waves + RSA.

test_that("the closed loop produces physiological pressure with emergent variability", {
  cr <- cardioRespiratory(duration = 300, hr0 = 70, resp_rate = 12, seed = 1)
  expect_s3_class(cr, "cardio_respiratory")
  expect_true(all(is.finite(c(cr$systolic, cr$diastolic, cr$rr))))
  expect_gt(mean(cr$systolic), 110); expect_lt(mean(cr$systolic), 135)
  expect_gt(mean(cr$diastolic), 70); expect_lt(mean(cr$diastolic), 92)
  # blood pressure carries a Mayer (LF) oscillation and a respiratory (HF) oscillation
  expect_gt(cr$bpv$lf_peak, 0.06); expect_lt(cr$bpv$lf_peak, 0.15)     # Mayer ~0.1 Hz
  expect_equal(cr$bpv$hf_peak, 12 / 60, tolerance = 0.03)              # Traube-Hering at breathing rate
  expect_equal(cr$hrv$hf_peak, 12 / 60, tolerance = 0.03)              # RSA at breathing rate
  # baroreflex sensitivity recovered from the beat series (~ g_vagal in ms/mmHg)
  expect_gt(cr$brs, 5); expect_lt(cr$brs, 20)
  expect_output(print(cr), "Mayer")
})

test_that("the Mayer wave is intrinsic to the baroreflex loop, not driven by breathing", {
  on  <- cardioRespiratory(duration = 300, resp_rate = 12, seed = 1)
  off <- cardioRespiratory(duration = 300, respiration = FALSE, seed = 1)
  # with respiration switched off the Mayer (LF) oscillation of pressure survives ...
  expect_gt(off$bpv$lf_power, 0.5 * on$bpv$lf_power)
  expect_gt(off$bpv$lf_peak, 0.06); expect_lt(off$bpv$lf_peak, 0.15)
  # ... but the respiratory-frequency heart-rate oscillation (RSA) collapses
  expect_lt(off$hrv$hf_power, 0.4 * on$hrv$hf_power)
})

test_that("slower breathing moves the respiratory (HF) peak down", {
  fast <- cardioRespiratory(duration = 300, resp_rate = 12, seed = 1)   # 0.20 Hz
  slow <- cardioRespiratory(duration = 300, resp_rate = 9,  seed = 1)   # 0.15 Hz
  expect_lt(slow$hrv$hf_peak, fast$hrv$hf_peak)
  expect_equal(slow$bpv$hf_peak, 9 / 60, tolerance = 0.03)
})

test_that("cardioRespiratory is reproducible under a fixed seed", {
  expect_identical(cardioRespiratory(duration = 120, seed = 7)$systolic,
                   cardioRespiratory(duration = 120, seed = 7)$systolic)
})
