# Forward sensor model: optical motion-capture markers (dropout / occlusion).

test_that("marker occlusion matches the expected gap fraction and is contiguous", {
  pos <- matrix(0, 5000, 2)                                 # clean = 0, so noise is recoverable
  m <- markerMeasure(pos, markerSensor(noise = 0.003, dropout_prob = 0.02, dropout_len = 10), seed = 3)
  expect_s3_class(m, "marker_measurement")
  p <- 0.02; L <- 10; expected <- p * L / (1 + p * L)       # stationary occluded fraction
  expect_lt(abs(m$gap_fraction - expected), 0.05)
  # occlusions are runs, not isolated frames (mean run length ~ dropout_len)
  runs <- rle(!m$visible); expect_gt(mean(runs$lengths[runs$values]), 4)
  # positional noise on visible frames recovers the sensor noise
  expect_equal(stats::sd((m$position - pos)[m$visible, ]), 0.003, tolerance = 0.1)
  expect_output(print(m), "occluded")
})

test_that("gap fraction increases with dropout probability", {
  pos <- matrix(0, 5000, 2)
  g1 <- markerMeasure(pos, markerSensor(dropout_prob = 0.01), seed = 1)$gap_fraction
  g2 <- markerMeasure(pos, markerSensor(dropout_prob = 0.05), seed = 1)$gap_fraction
  expect_gt(g2, g1)
})

test_that("fillMarkerGaps recovers the occluded trajectory", {
  set.seed(1); pos <- cbind(cumsum(rnorm(2000, 0, 0.01)), cumsum(rnorm(2000, 0, 0.01)))
  m <- markerMeasure(pos, markerSensor(noise = 0.001, dropout_prob = 0.03, dropout_len = 10), seed = 2)
  filled <- fillMarkerGaps(m)
  expect_false(any(is.na(filled)))                          # gaps are filled
  expect_lt(sqrt(mean((filled - pos)^2)), 0.03)             # close to the true trajectory
})

test_that("markerMeasure is reproducible under a fixed seed", {
  pos <- matrix(rnorm(600), 300, 2)
  a <- markerMeasure(pos, markerSensor(dropout_prob = 0.03), seed = 7)
  b <- markerMeasure(pos, markerSensor(dropout_prob = 0.03), seed = 7)
  expect_identical(a$visible, b$visible)
  expect_equal(a$position, b$position)
})
