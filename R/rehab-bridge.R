# Layer 5b -- close the clinical loop against the real-world rehabilitation
# evaluation framework (PhysioRehab).
#
# The twin PREDICTS an intervention effect; the clinic then collects a real
# single-case series (a baseline phase, then an intervention phase). This bridges
# the two: it hands that real series to PhysioRehab's single-case experimental
# design (SCED) engine to get the clinically-standard verdict -- the NAP / Tau
# effect size, whether the change exceeds the minimal clinically important
# difference (MCID), and the corresponding ICF category -- and then RECONCILES the
# twin's prediction with what actually happened. This is how an in-silico
# prediction is held to account against real, clinically-graded evidence.
#
# PhysioRehab (and its SCED engine PhysioAppKit) is an optional dependency
# (Suggests); this bridge is inert without it.

#' Evaluate a twin's intervention prediction against a real single-case series
#'
#' Connects the digital-twin decision loop to the PhysioRehab clinical-reasoning
#' framework. A real single-case series -- outcome measurements in a `baseline`
#' phase then an `intervention` phase -- is assembled into a `rehab_episode` and
#' analysed with PhysioRehab's single-case experimental design engine: the NAP and
#' Tau effect sizes, whether the change exceeds the `mcid`, and the ICF category of
#' the `measure`. If the twin's `predicted_delta` is supplied, it is reconciled
#' with the clinically-observed change (do they agree in sign, and on the MCID
#' verdict?), closing the predict -> deliver -> evaluate loop against real evidence.
#'
#' @param baseline Numeric outcome measurements in the baseline phase.
#' @param intervention Numeric outcome measurements in the intervention phase.
#' @param predicted_delta Optional: the twin's predicted change in the outcome
#'   (e.g. from [scoreInterventions()] / [evaluateIntervention()]), for reconciliation.
#' @param measure Outcome name (a PhysioRehab measure, e.g. `"gait_speed"`,
#'   `"cadence"`; drives the ICF mapping).
#' @param mcid Minimal clinically important difference (threshold, outcome units).
#' @param direction `"increase"` or `"decrease"` -- the improving direction.
#' @param condition,patient_id Episode metadata.
#' @param start_date First session date (`Date` or a string coercible by [as.Date()]).
#' @return a `rehab_evaluation`: the built `episode`, the PhysioRehab `sced` result,
#'   the `icf` code(s), and -- when `predicted_delta` is given -- a `reconciliation`
#'   (`predicted`, `observed`, `sign_agree`, `mcid_agree`, `error`).
#' @seealso [evaluateIntervention()], [clinicalDecision()]
#' @export
#' @examples
#' if (requireNamespace("PhysioRehab", quietly = TRUE)) {
#'   set.seed(1)
#'   base <- rnorm(5, 0.80, 0.03); post <- rnorm(6, 1.00, 0.03)
#'   rehabEvaluate(base, post, predicted_delta = 0.18, measure = "gait_speed", mcid = 0.16)
#' }
rehabEvaluate <- function(baseline, intervention, predicted_delta = NULL,
                          measure = "gait_speed", mcid = 0.16, direction = "increase",
                          condition = "in-silico twin", patient_id = "twin",
                          start_date = "2026-01-01") {
  if (!requireNamespace("PhysioRehab", quietly = TRUE))
    stop("rehabEvaluate() connects the twin to PhysioRehab's SCED/MCID/ICF workflow; ",
         "install the 'PhysioRehab' package (listed in Suggests).", call. = FALSE)
  baseline <- as.numeric(baseline); intervention <- as.numeric(intervention)
  named <- measure %in% c("gait_speed", "cadence", "symmetry")   # PhysioRehab named columns
  d0 <- as.Date(start_date); k <- 0L
  ep <- PhysioRehab::new_rehab_episode(patient_id = patient_id, condition = condition, design = "AB")
  add <- function(ep, val, phase) {
    k <<- k + 1L
    a <- list(ep = ep, date = d0 + k, phase = phase)
    if (named) a[[measure]] <- val else a$measures <- stats::setNames(list(val), measure)
    do.call(PhysioRehab::add_session, a)
  }
  for (v in baseline) ep <- add(ep, v, "baseline")
  for (v in intervention) ep <- add(ep, v, "intervention")
  s <- PhysioRehab::sced_analyze(ep, measure = measure, mcid = mcid, direction = direction)
  icf <- tryCatch(PhysioRehab::icf_for_measure(measure), error = function(e) character(0))
  recon <- NULL
  if (!is.null(predicted_delta)) {
    obs <- s$delta
    recon <- list(predicted = predicted_delta, observed = obs,
                  sign_agree = sign(predicted_delta) == sign(obs),
                  mcid_agree = (predicted_delta >= mcid) == isTRUE(s$exceeds_mcid),
                  error = obs - predicted_delta)
  }
  structure(list(episode = ep, sced = s, icf = icf, measure = measure, mcid = mcid,
                 predicted_delta = predicted_delta, reconciliation = recon),
            class = "rehab_evaluation")
}

#' @export
print.rehab_evaluation <- function(x, ...) {
  s <- x$sced
  cat(sprintf("Rehab evaluation -- %s (ICF %s)\n", x$measure,
              if (length(x$icf)) paste(x$icf, collapse = "/") else "unmapped"))
  cat(sprintf("  SCED: NAP %.2f, Tau %.2f, change %+.3f -> %s (MCID %.3f)\n",
              s$nap, s$tau, s$delta,
              if (isTRUE(s$exceeds_mcid)) "clinically meaningful (>= MCID)" else "below MCID", x$mcid))
  if (!is.null(x$reconciliation)) {
    r <- x$reconciliation
    cat(sprintf("  twin vs reality: predicted %+.3f, observed %+.3f -- sign %s, MCID call %s\n",
                r$predicted, r$observed, if (r$sign_agree) "agree" else "DISAGREE",
                if (r$mcid_agree) "agree" else "DISAGREE"))
  }
  invisible(x)
}
