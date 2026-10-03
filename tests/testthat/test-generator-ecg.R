# ECGSYN mechanistic cardiac generator.

# refractory R-peak detector (200 ms) for checking the generated rate
detect_r <- function(z, fs, refr = 0.2) {
  th <- 0.6 * max(z); n <- length(z); pk <- integer(0); last <- -Inf
  for (i in 2:(n - 1)) if (z[i] > th && z[i] >= z[i - 1] && z[i] > z[i + 1] &&
                           (i - last) > refr * fs) { pk <- c(pk, i); last <- i }
  pk
}

test_that("ECGSYN produces a physiological waveform at the requested rate", {
  e <- ecgsyn(duration = 12, hr = 70, sfecg = 256, hrv_sd = 0.03, seed = 1)
  expect_s3_class(e, "ecgsyn")
  expect_equal(range(e$ecg), c(-0.4, 1.2), tolerance = 1e-6)   # R up (1.2), Q/S down (-0.4)
  r <- detect_r(e$ecg, 256)
  expect_equal(length(r) / 12 * 60, 70, tolerance = 12)        # ~ 70 bpm
  expect_equal(mean(e$rr), 60 / 70, tolerance = 0.05)
})

test_that("morphology has P/QRS/T with R the tallest and a negative QS", {
  e <- ecgsyn(duration = 10, hr = 60, sfecg = 512, hrv_sd = 0, seed = 2)
  fs <- 512; r <- detect_r(e$ecg, fs)
  # average beat aligned on R
  win <- round(-0.25 * fs):round(0.35 * fs)
  beats <- do.call(rbind, lapply(r[r + max(win) <= length(e$ecg) & r + min(win) >= 1],
                                 function(i) e$ecg[i + win]))
  avg <- colMeans(beats); tR <- which(win == 0)
  expect_equal(max(avg), avg[tR], tolerance = 1e-6)            # R is the peak of the beat
  baseline <- avg[which.min(abs(win + 0.09 * fs))]            # PR segment (before Q)
  s_wave   <- min(avg[win >= 0 & win <= 0.08 * fs])           # S trough after R
  expect_lt(s_wave, baseline)                                 # Q/S dips below baseline
  st  <- avg[which.min(abs(win - 0.13 * fs))]                  # ST segment
  tpk <- max(avg[win >= 0.15 * fs & win <= 0.35 * fs])         # T-wave region
  expect_gt(tpk, st)                                           # T wave rises above the ST segment
})

test_that("HRV: fixed rate has near-zero SDNN; hrv_sd raises it", {
  fixed <- ecgsyn(duration = 20, hr = 60, hrv_sd = 0, seed = 3)
  vary  <- ecgsyn(duration = 20, hr = 60, hrv_sd = 0.06, seed = 3)
  expect_lt(sd(fixed$rr), 1e-6)
  expect_gt(sd(vary$rr), sd(fixed$rr))
  expect_gt(sd(vary$rr), 0.03)
  expect_output(print(vary), "ECGSYN")
})

test_that("heart rate maps to the RR interval", {
  slow <- ecgsyn(duration = 10, hr = 50, hrv_sd = 0)
  fast <- ecgsyn(duration = 10, hr = 100, hrv_sd = 0)
  expect_equal(mean(slow$rr), 60 / 50, tolerance = 0.02)
  expect_equal(mean(fast$rr), 60 / 100, tolerance = 0.02)
  expect_gt(length(fast$beats), length(slow$beats))
})
