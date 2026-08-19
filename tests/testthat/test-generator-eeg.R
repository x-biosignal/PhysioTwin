# Physiological EEG generator: the Jansen-Rit neural-mass model.

test_that("Jansen-Rit produces a self-organised alpha rhythm with default drive", {
  e <- jansenRit(duration = 8, seed = 1)
  expect_s3_class(e, "jansenRit")
  expect_equal(length(e$eeg), 8 * 256)
  expect_true(all(is.finite(e$eeg)))                     # no divergence
  expect_gt(sd(e$eeg), 0.2)                              # a real oscillation, not flat
  # the emergent rhythm is in the alpha band ...
  expect_gt(e$dominant_hz, 8); expect_lt(e$dominant_hz, 13)
  expect_identical(e$band, "alpha")
  # ... and alpha is the dominant band of the spectrum
  bp <- bandPower(e)
  expect_lte(sum(bp), 1 + 1e-9)            # relative distribution (sub-delta excluded)
  expect_gt(sum(bp), 0.9)                  # almost all power in the named bands
  expect_identical(names(which.max(bp)), "alpha")
  expect_output(print(e), "Jansen-Rit")
})

test_that("dropping the background drive collapses the alpha limit cycle to rest", {
  alpha <- jansenRit(duration = 8, seed = 1)
  rest  <- jansenRit(duration = 8, p_mean = 50, p_sd = 5, seed = 1)
  post  <- function(x) { s <- x$eeg[floor(length(x$eeg) / 3):length(x$eeg)]; sd(s) }
  # below the oscillation threshold the column settles to a low-amplitude fixed point
  expect_lt(post(rest), 0.2 * post(alpha))
  expect_gt(post(alpha), 0.3)
})

test_that("slower membrane time constants give a slower rhythm", {
  alpha <- jansenRit(duration = 8, seed = 1)                 # a = 100, b = 50
  slow  <- jansenRit(duration = 8, a = 70, b = 35, seed = 1) # slower synaptic dynamics
  expect_lt(slow$dominant_hz, alpha$dominant_hz - 2)         # a clear downward shift
  expect_lt(slow$dominant_hz, 8)                             # out of alpha, into theta/delta
})

test_that("jansenRit is reproducible under a fixed seed", {
  expect_identical(jansenRit(duration = 3, seed = 42)$eeg,
                   jansenRit(duration = 3, seed = 42)$eeg)
  # different seeds give different background-driven realisations
  expect_false(identical(jansenRit(duration = 3, seed = 1)$eeg,
                         jansenRit(duration = 3, seed = 2)$eeg))
})

test_that("aperiodicNoise has a 1/f^exponent spectrum and target SD", {
  x <- aperiodicNoise(4096, exponent = 1.5, sd = 2, seed = 1)
  expect_equal(sd(x), 2, tolerance = 0.05)
  expect_lt(abs(mean(x)), 0.05)
  sp <- spec.pgram(x, plot = FALSE, taper = 0); f <- sp$freq; ok <- f > 0.01 & f < 0.4
  slope <- coef(lm(log(sp$spec[ok]) ~ log(f[ok])))[2]
  expect_lt(abs(slope - (-1.5)), 0.3)                       # log-log slope ~ -exponent
  expect_identical(aperiodicNoise(256, seed = 7), aperiodicNoise(256, seed = 7))
})

test_that("a 1/f background makes the Jansen-Rit spectrum realistically broadband", {
  pure <- jansenRit(duration = 8, aperiodic = 0, seed = 1)
  real <- jansenRit(duration = 8, aperiodic = 1, seed = 1)
  # the pure column over-concentrates alpha; the 1/f background spreads it toward realistic
  expect_gt(bandPower(pure)[["alpha"]], 0.95)
  expect_lt(bandPower(real)[["alpha"]], 0.8)
  expect_gt(bandPower(real)[["alpha"]], 0.4)                # still the largest band, not overwhelming
  # the alpha peak survives the 1/f floor (the dominant-rhythm search is above 3 Hz)
  expect_identical(real$band, "alpha")
  expect_gt(real$dominant_hz, 8); expect_lt(real$dominant_hz, 13)
})
