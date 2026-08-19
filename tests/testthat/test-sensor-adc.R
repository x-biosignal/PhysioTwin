# ADC artefacts: quantisation, resampling, jitter.

test_that("quantisation error matches the theoretical LSB/sqrt(12)", {
  sig <- sin(seq(0, 40, length.out = 8000))
  qn <- quantizationNoise(sig, bits = 12)
  expect_lt(abs(qn[["empirical_sd"]] - qn[["theoretical_sd"]]) / qn[["theoretical_sd"]], 0.2)
  expect_equal(unname(qn[["theoretical_power"]]), unname(qn[["lsb"]])^2 / 12, tolerance = 1e-9)
  # more bits -> smaller quantisation error (monotone)
  expect_gt(quantizationNoise(sig, 8)[["empirical_sd"]], quantizationNoise(sig, 12)[["empirical_sd"]])
  expect_gt(quantizationNoise(sig, 12)[["empirical_sd"]], quantizationNoise(sig, 16)[["empirical_sd"]])
})

test_that("adcQuantize limits the number of levels and stays in range", {
  sig <- sin(seq(0, 6, length.out = 1000))
  q <- adcQuantize(sig, bits = 6, range = c(-1, 1))
  expect_lte(length(unique(round(q, 9))), 2^6)
  expect_true(all(q >= -1 & q <= 1))
})

test_that("adcSample resamples to the new rate and preserves a low-frequency signal", {
  sig <- sin(2 * pi * 2 * seq(0, 4, length.out = 4000))     # 2 Hz over 4 s at 1000 Hz
  d <- adcSample(sig, fs_in = 1000, fs_out = 250)
  expect_equal(d$fs, 250)
  expect_lt(abs(length(d$signal) - 1000), 5)                # ~ 4 s x 250 Hz
  expect_gt(cor(d$signal, sin(2 * pi * 2 * d$time)), 0.99)  # the 2 Hz content survives
})

test_that("dithered quantisation is reproducible under a fixed seed", {
  sig <- sin(seq(0, 10, length.out = 500))
  set.seed(1); a <- adcQuantize(sig, bits = 8, dither = TRUE)
  set.seed(1); b <- adcQuantize(sig, bits = 8, dither = TRUE)
  expect_identical(a, b)
})
