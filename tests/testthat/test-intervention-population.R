# Population-level in-silico intervention (an in-silico trial).

test_that("the in-silico trial scales with the intervention dose", {
  skip_on_cran()
  big   <- populationIntervention(n = 100, intervention = list(scale = list(strength = 1.6)), seed = 1)
  small <- populationIntervention(n = 100, intervention = list(scale = list(strength = 1.1)), seed = 1)
  expect_s3_class(big, "population_intervention")
  expect_equal(nrow(big$subjects), 100)
  expect_gt(big$effect["mean"], small$effect["mean"])       # bigger dose, bigger effect
  expect_gt(big$effect["mean"], 0)                          # strengthening increases ROM
  expect_gte(big$responder_rate, small$responder_rate)
  expect_gte(big$power, small$power)
  expect_lte(big$n_for_80, small$n_for_80)                  # fewer subjects needed
  expect_output(print(big), "in-silico trial")
})

test_that("the analytic paired-t power matches a Monte-Carlo trial", {
  skip_on_cran()
  pop <- populationIntervention(n = 25, intervention = list(scale = list(strength = 1.02)), seed = 2)
  dz <- pop$effect["cohen_dz"]
  set.seed(11)
  mc <- mean(replicate(3000, stats::t.test(stats::rnorm(25, dz, 1))$p.value < 0.05))
  expect_lt(abs(pop$power - mc), 0.05)
})

test_that("power increases with sample size (noncentral-t)", {
  skip_on_cran()
  pop <- populationIntervention(n = 100, intervention = list(scale = list(strength = 1.02)), seed = 1)
  dz <- pop$effect["cohen_dz"]
  pw <- function(nn) { tc <- qt(0.975, nn - 1); ncp <- dz * sqrt(nn)
    pt(-tc, nn - 1, ncp) + pt(tc, nn - 1, ncp, lower.tail = FALSE) }
  expect_lt(pw(10), pw(100))                                # power rises with n
  expect_equal(pop$power, pw(100), tolerance = 1e-6)        # reported power uses the formula
  expect_true(is.na(pop$n_for_80) || pop$n_for_80 >= 2)
})

test_that("a custom outcome function works and results are reproducible", {
  skip_on_cran()
  peakacc <- function(tr) max(abs(tr$alpha))
  a <- populationIntervention(n = 40, outcome = peakacc,
                              intervention = list(scale = list(strength = 1.3)), seed = 5)
  b <- populationIntervention(n = 40, outcome = peakacc,
                              intervention = list(scale = list(strength = 1.3)), seed = 5)
  expect_equal(a$subjects, b$subjects)
  expect_identical(a$outcome, "custom")
})
