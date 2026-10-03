# Layer 4b -- population-level in-silico intervention (an in-silico trial).
#
# A single-subject in-silico intervention answers "how would THIS twin change?".
# A trial asks a population question: sample a cohort of twins from a parameter
# distribution, apply the same intervention to each, and characterise the
# distribution of responses -- the mean effect and its size, the responder rate
# against a clinically important difference, and the statistical power (and the
# sample size for 80% power) of a paired trial that would detect it. This turns the
# twin into an in-silico trial simulator for designing real studies.

# paired-t power from the (paired) effect size d_z at sample size nn
.paired_power <- function(dz, nn, alpha = 0.05) {
  if (nn < 2) return(0)
  tc <- stats::qt(1 - alpha / 2, nn - 1); ncp <- dz * sqrt(nn)
  stats::pt(-tc, nn - 1, ncp) + stats::pt(tc, nn - 1, ncp, lower.tail = FALSE)
}
# smallest n giving >= target power (binary search)
.n_for_power <- function(dz, target = 0.8, nmax = 1e5L) {
  if (.paired_power(dz, nmax) < target) return(NA_integer_)
  lo <- 2L; hi <- nmax
  while (lo < hi) { mid <- (lo + hi) %/% 2L
    if (.paired_power(dz, mid) >= target) hi <- mid else lo <- mid + 1L }
  lo
}

#' In-silico intervention trial over a population of twins
#'
#' Samples a cohort of movement twins from a parameter distribution, applies the
#' same intervention to each (via [insilicoIntervention()]), and summarises the
#' response as a paired in-silico trial: the mean effect and its paired effect size
#' (Cohen's `d_z`), the responder rate against a minimal clinically important
#' difference, the statistical power at this sample size and the sample size needed
#' for 80% power.
#'
#' @param n Number of subjects (twins) in the cohort.
#' @param params Named list of `c(mean, sd)` population distributions for
#'   [limbTwin()] parameters (sampled Gaussian, truncated positive).
#' @param intervention List with `scale` and/or `set` (as [insilicoIntervention()]).
#' @param outcome Outcome to track: one of `"rom"`, `"peak_velocity"`,
#'   `"mean_speed"`, `"peak_accel"`, or a `function(trajectory)` returning a scalar.
#' @param mcid Minimal clinically important difference (responder threshold);
#'   defaults to 0.1 x the baseline outcome SD.
#' @param duration,dt,theta0,omega0 Simulation settings (as [simulateTwin()]).
#' @param seed Optional RNG seed.
#' @return a `population_intervention`: `subjects` (per-subject params, `pre`,
#'   `post`, `delta`), `effect` (mean delta, SD, Cohen's `d_z`), `responder_rate`,
#'   `mcid`, `power` and `n_for_80`.
#' @seealso [insilicoIntervention()], [multimodalTwin()], [generateTrainingData()]
#' @export
#' @examples
#' \donttest{
#' pop <- populationIntervention(n = 100, intervention = list(scale = list(strength = 1.4)),
#'                               seed = 1)
#' pop$effect
#' }
populationIntervention <- function(n = 200,
                                   params = list(strength = c(5, 1), damping = c(0.7, 0.2),
                                                 stiffness = c(2.5, 0.6)),
                                   intervention = list(scale = list(strength = 1.4)),
                                   outcome = "rom", mcid = NULL,
                                   duration = 6, dt = 0.005, theta0 = 0.1, omega0 = 0,
                                   seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  ofun <- if (is.function(outcome)) outcome else {
    outcome <- match.arg(outcome, c("rom", "peak_velocity", "mean_speed", "peak_accel"))
    function(tr) .twin_outcomes(tr)[[outcome]]
  }
  pnames <- names(params)
  draws <- vapply(params, function(ms) pmax(1e-3, stats::rnorm(n, ms[1], ms[2])), numeric(n))
  draws <- matrix(draws, n, length(pnames), dimnames = list(NULL, pnames))
  scale <- intervention$scale %||% list(); set <- intervention$set %||% list()
  pre <- post <- numeric(n)
  for (i in seq_len(n)) {
    tw <- do.call(limbTwin, as.list(draws[i, ]))
    iv <- insilicoIntervention(tw, scale = scale, set = set, duration = duration,
                               dt = dt, theta0 = theta0, omega0 = omega0)
    pre[i] <- ofun(iv$before); post[i] <- ofun(iv$after)
  }
  delta <- post - pre
  if (is.null(mcid)) mcid <- 0.1 * stats::sd(pre)
  sdd <- stats::sd(delta)                                  # d_z is 0 for a true null, +/-Inf
  dz <- if (sdd < 1e-12) if (abs(mean(delta)) < 1e-12) 0 else sign(mean(delta)) * Inf else mean(delta) / sdd
  subjects <- data.frame(draws, pre = pre, post = post, delta = delta)
  structure(list(subjects = subjects,
                 effect = c(mean = mean(delta), sd = stats::sd(delta), cohen_dz = dz),
                 responder_rate = mean(delta > mcid), mcid = mcid,
                 power = .paired_power(dz, n), n_for_80 = .n_for_power(dz),
                 outcome = if (is.character(outcome)) outcome else "custom", n = n),
            class = "population_intervention")
}

#' @export
print.population_intervention <- function(x, ...) {
  cat(sprintf("Population in-silico trial -- n=%d, outcome '%s'\n", x$n, x$outcome))
  cat(sprintf("  mean effect %.3f (d_z %.2f), responders %.0f%% (MCID %.3f)\n",
              x$effect["mean"], x$effect["cohen_dz"], 100 * x$responder_rate, x$mcid))
  cat(sprintf("  power at n=%d: %.2f;  n for 80%% power: %s\n",
              x$n, x$power, ifelse(is.na(x$n_for_80), ">1e5", x$n_for_80)))
  invisible(x)
}

# NULL-coalescing helper (base R has none)
`%||%` <- function(a, b) if (is.null(a)) b else a
