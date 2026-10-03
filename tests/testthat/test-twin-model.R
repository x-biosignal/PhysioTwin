# Layer 1: the mechanistic movement twin.

test_that("simulateTwin integrates a bounded rhythmic movement", {
  tw <- limbTwin(strength = 5, damping = 0.4, stiffness = 3, freq = 1)
  sim <- simulateTwin(tw, duration = 6, dt = 0.005)
  expect_s3_class(sim, "twin_trajectory")
  expect_equal(length(sim$theta), length(sim$time))
  expect_true(all(is.finite(sim$theta)))
  expect_lt(diff(range(sim$theta)), 6)                 # bounded (not diverging)
  expect_gt(diff(range(sim$theta)), 0.05)              # actually moves
})

test_that("more muscle strength increases range of motion", {
  weak   <- simulateTwin(limbTwin(strength = 3, damping = 0.4), duration = 8, dt = 0.005)
  strong <- simulateTwin(limbTwin(strength = 9, damping = 0.4), duration = 8, dt = 0.005)
  expect_gt(diff(range(strong$theta)), diff(range(weak$theta)))
})

test_that("more damping reduces peak angular velocity", {
  lo <- simulateTwin(limbTwin(strength = 6, damping = 0.2), duration = 8, dt = 0.005)
  hi <- simulateTwin(limbTwin(strength = 6, damping = 1.2), duration = 8, dt = 0.005)
  expect_lt(max(abs(hi$omega)), max(abs(lo$omega)))
})

test_that("horizontal-plane (g = 0) movement runs and prints", {
  sim <- simulateTwin(limbTwin(g = 0, strength = 4), duration = 3, dt = 0.01)
  expect_true(all(is.finite(sim$omega)))
  expect_output(print(sim), "Twin trajectory")
  expect_output(print(limbTwin()), "Movement twin")
})

test_that("inertia must be positive", {
  expect_error(limbTwin(inertia = 0))
})
