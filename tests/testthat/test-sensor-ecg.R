# Forward sensor model: ECG lead / electrode.

test_that("the recorded lead carries baseline-wander and powerline peaks", {
  e <- ecgsyn(duration = 20, seed = 1)
  rec <- ecgMeasure(e$ecg, ecgLead(baseline_amp = 0.2, baseline_hz = 0.3,
                                   powerline_hz = 50, powerline_amp = 0.05, fs = e$sfecg), seed = 1)
  expect_s3_class(rec, "ecg_measurement")
  n <- length(rec$signal); f <- (0:(n - 1)) * e$sfecg / n; p <- Mod(stats::fft(rec$signal))^2
  # baseline-wander peak near 0.3 Hz
  lo <- f > 0.1 & f < 1; expect_lt(abs(f[lo][which.max(p[lo])] - 0.3), 0.1)
  # powerline peak near 50 Hz
  b50 <- which.min(abs(f - 50)); expect_gt(p[b50], 5 * stats::median(p[f > 40 & f < 60]))
  expect_output(print(rec), "ECG measurement")
})

test_that("R-peaks survive baseline wander and contamination", {
  e <- ecgsyn(duration = 20, seed = 1)
  rec <- ecgMeasure(e$ecg, ecgLead(baseline_amp = 0.2, baseline_hz = 0.3,
                                   powerline_amp = 0.05, noise = 0.02, fs = e$sfecg), seed = 1)
  rc <- detectRpeaks(e$ecg, e$sfecg); rr <- detectRpeaks(rec$signal, e$sfecg)
  expect_equal(length(rr), length(rc))                     # same number of beats recovered
  expect_lt(stats::median(vapply(rc, function(k) min(abs(rr - k)), numeric(1))), 3)  # aligned
})

test_that("SNR degrades monotonically as muscle noise rises", {
  e <- ecgsyn(duration = 20, seed = 1)
  snr <- vapply(c(0.01, 0.05, 0.15), function(nz) {
    r <- ecgMeasure(e$ecg, ecgLead(noise = nz, baseline_amp = 0.1, fs = e$sfecg), seed = 2)
    10 * log10(stats::var(e$ecg) / stats::var(r$signal - e$ecg))
  }, numeric(1))
  expect_true(all(diff(snr) < 0))
})

test_that("ecgMeasure is reproducible under a fixed seed", {
  e <- ecgsyn(duration = 6, seed = 1)
  expect_equal(ecgMeasure(e$ecg, ecgLead(), seed = 4)$signal,
               ecgMeasure(e$ecg, ecgLead(), seed = 4)$signal)
})
