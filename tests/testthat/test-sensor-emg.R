# Forward sensor model: surface-EMG electrode / amplifier.

emg_spectrum <- function(x, fs) {
  n <- length(x); list(f = (0:(n - 1)) * fs / n, p = Mod(stats::fft(x))^2)
}

test_that("the raw sEMG carries a powerline peak at the mains frequency", {
  act <- pmax(0, sin(seq(0, 8 * pi, length.out = 4000)))^1.5
  emg <- emgMeasure(act, emgElectrode(fs = 1000, powerline_hz = 50, powerline_amp = 0.05), seed = 1)
  expect_s3_class(emg, "emg_measurement")
  s <- emg_spectrum(emg$raw, 1000); b50 <- which.min(abs(s$f - 50))
  nbr <- stats::median(s$p[c((b50 - 20):(b50 - 5), (b50 + 5):(b50 + 20))])
  expect_gt(s$p[b50], 5 * nbr)                              # sharp mains peak
  expect_output(print(emg), "EMG measurement")
})

test_that("the EMG envelope recovers the activation", {
  act <- pmax(0, sin(seq(0, 8 * pi, length.out = 4000)))^1.5
  emg <- emgMeasure(act, emgElectrode(fs = 1000), seed = 1)
  expect_gt(cor(emgEnvelope(emg$clean, 1000), act), 0.8)    # amplitude tracks the drive
})

test_that("the clean EMG is band-limited to the electrode passband", {
  act <- rep(1, 4000)                                       # constant drive -> broadband source
  emg <- emgMeasure(act, emgElectrode(fs = 1000, bandwidth = c(20, 450)), seed = 1)
  s <- emg_spectrum(emg$clean, 1000)
  ff <- pmin(s$f, 1000 - s$f)                               # fold the two-sided spectrum to [0, fs/2]
  inband <- ff >= 20 & ff <= 450
  expect_gt(sum(s$p[inband]) / sum(s$p), 0.9)
})

test_that("emgMeasure is reproducible under a fixed seed", {
  act <- pmax(0, sin(seq(0, 6, length.out = 1000)))
  expect_equal(emgMeasure(act, emgElectrode(), seed = 9)$raw,
               emgMeasure(act, emgElectrode(), seed = 9)$raw)
})
