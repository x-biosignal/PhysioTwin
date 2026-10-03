# Cardiovascular generators: Windkessel arterial pressure + baroreflex.

test_that("Windkessel produces a physiological pressure waveform", {
  bp <- windkessel(duration = 12, hr = 70)
  expect_s3_class(bp, "windkessel")
  expect_equal(length(bp$pressure), 12 * 250)
  expect_true(all(is.finite(bp$pressure)))
  # normal adult ranges
  expect_gt(bp$systolic, 110); expect_lt(bp$systolic, 140)
  expect_gt(bp$diastolic, 65); expect_lt(bp$diastolic, 90)
  expect_gt(bp$systolic, bp$diastolic)
  expect_equal(bp$pp, bp$systolic - bp$diastolic)
  # mean arterial pressure tracks cardiac output x resistance
  expect_equal(bp$map, (80 * 70 / 60) * 1.0, tolerance = 0.1)
  expect_output(print(bp), "Windkessel")
})

test_that("the diastolic run-off decays with time constant RC", {
  bp <- windkessel(duration = 14, hr = 60, R = 1.0, C = 1.9)   # RC = 1.9 s
  i <- bp$time > 12 & bp$time < 13                             # one settled beat
  p <- bp$pressure[i]; tt <- (seq_along(p) - 1) / 250
  dia <- tt > 0.3 + 0.05                                       # pure diastole (no inflow)
  tau <- as.numeric(-1 / coef(lm(log(p[dia]) ~ tt[dia]))[2])
  expect_equal(tau, bp$tau, tolerance = 0.15)                  # decays toward 0 with RC
})

test_that("resistance raises mean pressure and low compliance widens pulse pressure", {
  base  <- windkessel(duration = 12)
  hyper <- windkessel(duration = 12, R = 1.4)                  # vasoconstriction
  stiff <- windkessel(duration = 12, C = 0.9)                  # arterial stiffening
  expect_gt(hyper$map, base$map + 10)                          # hypertension
  expect_gt(stiff$pp, base$pp + 10)                            # wider pulse pressure
})

test_that("the baroreflex buffers a pressure perturbation", {
  on  <- baroreflex(duration = 80, perturb_time = 40, perturb_R = 0.7, reflex = TRUE)
  off <- baroreflex(duration = 80, perturb_time = 40, perturb_R = 0.7, reflex = FALSE)
  expect_s3_class(on, "baroreflex")
  # with the reflex the post-perturbation mean-pressure error is much smaller
  expect_lt(on$regulated, off$regulated)
  expect_lt(on$regulated, 0.6 * off$regulated)
  # a fall in resistance (vasodilation) triggers a compensatory rise in heart rate
  base_hr <- on$hr[which.min(abs(on$beat_time - 38))]
  post_hr <- mean(tail(on$hr, 5))
  expect_gt(post_hr, base_hr)
  expect_output(print(on), "closed loop")
})

test_that("baroreflex sensitivity scales the strength of regulation", {
  weak   <- baroreflex(duration = 80, perturb_time = 40, gain = 0.005)
  strong <- baroreflex(duration = 80, perturb_time = 40, gain = 0.04)
  expect_lt(strong$regulated, weak$regulated)          # tighter reflex -> smaller error
})
