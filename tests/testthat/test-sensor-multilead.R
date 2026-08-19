# Multi-lead ECG from one cardiac source.

test_that("multi-lead ECG projects one source into correlated but distinct leads", {
  e <- ecgsyn(duration = 15, seed = 1)
  ml <- ecgMeasureLeads(e$ecg, ecgLeadSet(leads = c("II", "V2", "V5"), fs = e$sfecg), seed = 1)
  expect_s3_class(ml, "ecg_multilead")
  expect_equal(dim(ml$signals), c(length(e$ecg), 3))
  # every lead is strongly correlated with the cardiac source ...
  cs <- apply(ml$signals, 2, function(s) cor(s, e$ecg))
  expect_true(all(cs > 0.7))
  # ... but the leads are distinct (not identical amplitude/shape)
  pw <- cor(ml$signals)[upper.tri(diag(3))]
  expect_true(all(pw < 0.99))
  expect_output(print(ml), "Multi-lead")
})

test_that("R-peaks align across leads", {
  e <- ecgsyn(duration = 20, seed = 2)
  ml <- ecgMeasureLeads(e$ecg, ecgLeadSet(leads = c("II", "V2", "V5"), fs = e$sfecg), seed = 2)
  rp <- lapply(1:3, function(l) detectRpeaks(ml$signals[, l], e$sfecg))
  k <- min(vapply(rp, length, integer(1)))
  expect_gt(k, 10)                                       # a plausible number of beats detected
  expect_lt(max(abs(rp[[1]][seq_len(k)] - rp[[3]][seq_len(k)])), 5)   # II vs V5 within a few samples
})

test_that("ecgMeasureLeads is reproducible", {
  e <- ecgsyn(duration = 8, seed = 1); s <- ecgLeadSet(fs = e$sfecg)
  expect_identical(ecgMeasureLeads(e$ecg, s, seed = 3)$signals,
                   ecgMeasureLeads(e$ecg, s, seed = 3)$signals)
})
