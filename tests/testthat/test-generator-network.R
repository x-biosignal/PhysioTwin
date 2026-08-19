# Coordination / network generators: Kuramoto + Matsuoka CPG.

test_that("Kuramoto shows a synchronisation transition with coupling", {
  incoh <- kuramoto(n = 20, K = 0, duration = 20, dt = 0.02, seed = 1)
  part  <- kuramoto(n = 20, K = 1, duration = 20, dt = 0.02, seed = 1)
  sync  <- kuramoto(n = 20, K = 4, duration = 20, dt = 0.02, seed = 1)
  expect_s3_class(sync, "kuramoto")
  expect_true(all(c(incoh$mean_order, sync$mean_order) >= 0 &
                  c(incoh$mean_order, sync$mean_order) <= 1))
  expect_lt(incoh$mean_order, 0.45)                    # incoherent
  expect_gt(sync$mean_order, 0.8)                      # synchronised
  expect_gt(sync$mean_order, part$mean_order)          # more coupling -> more order
  expect_gt(part$mean_order, incoh$mean_order)
  expect_output(print(sync), "Kuramoto")
})

test_that("Kuramoto with identical frequencies fully synchronises", {
  k <- kuramoto(n = 12, K = 2, omega = rep(1, 12), duration = 15, dt = 0.02, seed = 2)
  expect_gt(k$mean_order, 0.95)
  expect_equal(dim(k$signals), c(750, 12))
})

test_that("Matsuoka CPG produces a sustained anti-phase rhythm", {
  cpg <- cpgMatsuoka(duration = 20, dt = 0.005)
  expect_s3_class(cpg, "cpg")
  expect_gt(diff(range(cpg$output)), 0.2)             # a real oscillation, not flat
  expect_gt(cpg$frequency, 0.05); expect_lt(cpg$frequency, 5)   # a plausible rhythm
  expect_lt(cor(cpg$y1, cpg$y2), -0.5)                # the two half-centres are anti-phase
  # the rhythm is sustained (late-window amplitude comparable to mid-window)
  n <- length(cpg$output)
  expect_gt(sd(cpg$output[round(0.8 * n):n]), 0.05)
  expect_output(print(cpg), "Matsuoka CPG")
})

test_that("faster CPG time constants give a faster rhythm", {
  slow <- cpgMatsuoka(tau = 0.5, tau2 = 1.0, duration = 25, dt = 0.005)
  fast <- cpgMatsuoka(tau = 0.15, tau2 = 0.3, duration = 25, dt = 0.005)
  expect_gt(fast$frequency, slow$frequency)
})
