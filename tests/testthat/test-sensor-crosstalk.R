# Multi-channel EMG array: crosstalk + amplifier chain.

test_that("crosstalk mixes channels and the matrix is well-formed", {
  expect_equal(crosstalkMatrix(4, 0), diag(4))
  X <- crosstalkMatrix(5, 0.2)
  expect_equal(rowSums(X), rep(1, 5))                       # energy-preserving
  expect_true(all(diag(X) > 0.5))                           # diagonal dominant
  set.seed(1)
  act <- matrix(pmax(0, sin(outer(seq(0, 10, length.out = 2500), c(1, 1.3, 1.7, 2.1)))), ncol = 4)
  off <- function(M) { m <- cor(M); mean(abs(m[upper.tri(m)])) }
  none <- emgArrayMeasure(act, emgArray(channels = 4, crosstalk = 0.0), seed = 1)
  lots <- emgArrayMeasure(act, emgArray(channels = 4, crosstalk = 0.3), seed = 1)
  expect_gt(off(lots$raw), off(none$raw) + 0.1)            # crosstalk raises cross-channel correlation
  expect_equal(dim(lots$raw), c(2500, 4))
  expect_output(print(lots), "EMG array")
})

test_that("a higher CMRR rejects more of the common-mode mains", {
  zero <- matrix(0, 4000, 4)                               # no muscle signal -> the 50 Hz peak IS the mains
  pk50 <- function(mm) { sp <- spec.pgram(mm$raw[, 1], plot = FALSE, taper = 0)
    sp$spec[which.min(abs(sp$freq * mm$fs - 50))] }
  lo <- emgArrayMeasure(zero, emgArray(cmrr_db = 40, powerline_amp = 0.5, noise = 1e-4), seed = 1)
  hi <- emgArrayMeasure(zero, emgArray(cmrr_db = 100, powerline_amp = 0.5, noise = 1e-4), seed = 1)
  expect_gt(pk50(lo), 100 * pk50(hi))                     # 60 dB more CMRR -> vastly smaller 50 Hz peak
})

test_that("the DC-blocking high-pass removes an electrode offset and is reproducible", {
  zero <- matrix(0, 3000, 3); s <- emgArray(channels = 3, offset = 2, noise = 0, powerline_amp = 0)
  m <- emgArrayMeasure(zero, s, seed = 1)
  expect_lt(abs(mean(m$raw[500:3000, 1])), 0.1)           # constant electrode offset removed
  expect_identical(m$raw, emgArrayMeasure(zero, s, seed = 1)$raw)
})
