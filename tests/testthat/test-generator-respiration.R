# Respiration + respiratory sinus arrhythmia (RSA) generators.

test_that("respiration produces an asymmetric breathing waveform at the set rate", {
  rp <- respiration(duration = 40, rate = 15, ie_ratio = 0.5)
  expect_s3_class(rp, "respiration")
  expect_equal(length(rp$volume), 40 * 100)
  expect_equal(mean(rp$volume), 0, tolerance = 1e-8)          # zero-mean
  # dominant frequency equals the breathing rate
  nc <- sum(abs(diff(sign(rp$volume))) > 0) / 2
  expect_equal(nc / max(rp$time), 15 / 60, tolerance = 0.03)
  # expiration-dominant: more samples falling than rising (I:E = 1:2)
  expect_lt(sum(diff(rp$volume) > 0), sum(diff(rp$volume) < 0))
  expect_output(print(rp), "Respiration")
})

test_that("RSA tachogram carries an HF peak at the respiratory rate", {
  h <- rsaTachogram(duration = 300, hr0 = 70, resp_rate = 15, seed = 1)
  expect_s3_class(h, "rsaTachogram")
  expect_equal(mean(h$hr), 70, tolerance = 2)                 # mean rate preserved
  # the high-frequency HRV peak sits at the breathing frequency (0.25 Hz)
  expect_equal(h$hf_peak, 15 / 60, tolerance = 0.04)
  expect_gt(h$hf_power, 0)
})

test_that("deeper respiratory modulation increases HF power", {
  shallow <- rsaTachogram(duration = 300, rsa = 0.06, seed = 1)
  deep    <- rsaTachogram(duration = 300, rsa = 0.20, seed = 1)
  expect_gt(deep$hf_power, shallow$hf_power)
})

test_that("slower breathing moves the HF peak to a lower frequency", {
  fast <- rsaTachogram(duration = 300, resp_rate = 15, seed = 1)   # 0.25 Hz
  slow <- rsaTachogram(duration = 300, resp_rate = 9,  seed = 1)   # 0.15 Hz
  expect_lt(slow$hf_peak, fast$hf_peak)
  expect_equal(slow$hf_peak, 9 / 60, tolerance = 0.04)
})

test_that("rsaTachogram is reproducible under a fixed seed", {
  expect_identical(rsaTachogram(duration = 120, seed = 42)$rr,
                   rsaTachogram(duration = 120, seed = 42)$rr)
})
