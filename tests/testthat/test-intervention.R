# Layer 4a: in-silico intervention.

test_that("strengthening increases range of motion, damping reduces peak velocity", {
  tw <- limbTwin(strength = 5, damping = 0.6, stiffness = 3)
  strengthen <- insilicoIntervention(tw, scale = list(strength = 1.4), duration = 8)
  expect_s3_class(strengthen, "insilico_intervention")
  rom <- strengthen$outcomes[strengthen$outcomes$outcome == "rom", ]
  expect_gt(rom$after, rom$before)                     # more ROM after strengthening
  expect_gt(rom$pct_change, 0)

  brace <- insilicoIntervention(tw, scale = list(damping = 2), duration = 8)
  pv <- brace$outcomes[brace$outcomes$outcome == "peak_velocity", ]
  expect_lt(pv$after, pv$before)                       # damping slows the movement
})

test_that("set overrides an absolute parameter value", {
  tw <- limbTwin(strength = 5)
  iv <- insilicoIntervention(tw, set = list(strength = 10), duration = 5)
  expect_equal(iv$twin$strength, 10)
  expect_output(print(iv), "In-silico intervention")
})

test_that("no change yields ~zero effect", {
  tw <- limbTwin()
  iv <- insilicoIntervention(tw, duration = 5)
  expect_equal(iv$outcomes$pct_change, rep(0, nrow(iv$outcomes)), tolerance = 1e-6)
})
