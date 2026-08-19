# Layer 5 -- clinical decision loop: personalise -> gate -> score -> reconcile ->
# evaluate.
#
# A clinician-facing loop built on the twin, with the honest-scope discipline that
# makes it safe: a personalised twin is only allowed to advise once it PASSES A
# VALIDATION GATE, every in-silico intervention is scored WITH its predictive
# uncertainty (never a point estimate), the mechanistic prediction is RECONCILED
# with external evidence rather than overriding it, and a real follow-up is fed
# back to check the prediction held. This is research tooling, not a medical
# device: predictions are conditional on the twin being validated for the patient.
#
# Pieces:
#   * intervention() / defaultInterventions() -- candidate interventions as named
#     parameter changes with burden/risk (and optional evidence) metadata.
#   * twinReadiness()   -- the validation gate: fit quality + parameter
#     identifiability (Fisher covariance) -> a verdict and a [0,1] confidence.
#   * scoreInterventions() -- run each candidate over the parameter posterior,
#     giving a predicted effect, a predictive interval, the probability of a
#     clinically important benefit, and a multi-criteria score (efficacy, certainty,
#     risk, burden) scaled by the twin's confidence.
#   * reconcileEvidence() -- precision-weighted combination of the mechanistic and
#     evidence-based effects, flagging disagreement.
#   * evaluateIntervention() -- did the observed follow-up fall in the predicted
#     interval? (closes the loop).
#   * clinicalDecision()  -- the whole loop end-to-end.
# Dependency-free base R; reuses .rmvn() and `%||%`.

# numeric Hessian by central finite differences
.num_hessian <- function(f, x, h = 1e-4) {
  n <- length(x); H <- matrix(0, n, n); f0 <- f(x)
  step <- function(i, s) { z <- x; z[i] <- z[i] + s; z }
  hs <- h * pmax(1, abs(x))
  for (i in seq_len(n)) for (j in i:n) {
    if (i == j) {
      H[i, i] <- (f(step(i, hs[i])) - 2 * f0 + f(step(i, -hs[i]))) / hs[i]^2
    } else {
      zpp <- x; zpp[i] <- zpp[i] + hs[i]; zpp[j] <- zpp[j] + hs[j]
      zpm <- x; zpm[i] <- zpm[i] + hs[i]; zpm[j] <- zpm[j] - hs[j]
      zmp <- x; zmp[i] <- zmp[i] - hs[i]; zmp[j] <- zmp[j] + hs[j]
      zmm <- x; zmm[i] <- zmm[i] - hs[i]; zmm[j] <- zmm[j] - hs[j]
      H[i, j] <- H[j, i] <- (f(zpp) - f(zpm) - f(zmp) + f(zmm)) / (4 * hs[i] * hs[j])
    }
  }
  H
}

# draw a parameter posterior (log-space) from a readiness gate; NULL cov -> point estimate
.posterior_draw <- function(readiness, n_mc) {
  if (!is.null(readiness) && all(is.finite(readiness$cov)))
    list(S = .rmvn(n_mc, log(readiness$estimates), readiness$cov), estnm = readiness$estimated)
  else list(S = NULL, estnm = character(0))
}

# outcome-change samples for one candidate over a shared parameter draw S (paired
# across candidates so their comparison is not confounded by sampling noise)
.delta_samples <- function(twin, cand, S, estnm, outcome, duration, dt, th0) {
  outfun <- function(tr) .twin_outcomes(tr)[[outcome]]
  nmc <- if (is.null(S)) 1L else nrow(S)
  d <- numeric(nmc)
  for (m in seq_len(nmc)) {
    tw <- twin; if (!is.null(S)) tw[estnm] <- as.list(exp(S[m, ]))
    iv <- insilicoIntervention(tw, scale = cand$scale, set = cand$set,
                               duration = duration, dt = dt, theta0 = th0, omega0 = 0)
    d[m] <- outfun(iv$after) - outfun(iv$before)
  }
  d
}

#' Define a candidate in-silico intervention
#'
#' A candidate intervention for the clinical loop: a named set of parameter
#' changes (multiplicative `scale` and/or absolute `set`, as [insilicoIntervention()])
#' plus its clinical `burden` and `risk`, and optionally an evidence-based effect.
#'
#' @param name Label.
#' @param scale,set Parameter changes (see [insilicoIntervention()]).
#' @param burden,risk Clinical burden and risk, each on a 0..1 scale.
#' @param evidence Optional `list(effect, se)`: an external (e.g. literature)
#'   estimate of the outcome change and its standard error, for reconciliation.
#' @return a `twin_intervention`.
#' @seealso [defaultInterventions()], [scoreInterventions()]
#' @export
#' @examples
#' intervention("strengthening", scale = list(strength = 1.4), burden = 0.4, risk = 0.1)
intervention <- function(name, scale = list(), set = list(), burden = 0, risk = 0,
                         evidence = NULL) {
  structure(list(name = name, scale = scale, set = set, burden = burden,
                 risk = risk, evidence = evidence), class = "twin_intervention")
}

#' @export
print.twin_intervention <- function(x, ...) {
  ch <- paste(c(sprintf("%s x%.2f", names(x$scale), unlist(x$scale)),
                sprintf("%s=%.2f", names(x$set), unlist(x$set))), collapse = ", ")
  cat(sprintf("Intervention '%s' -- %s (burden %.1f, risk %.1f)\n", x$name,
              if (nzchar(ch)) ch else "no change", x$burden, x$risk)); invisible(x)
}

#' A default library of movement-twin interventions
#'
#' A small library of candidate interventions for the single-joint movement twin:
#' strengthening, tone reduction, stretching, an orthosis, and a no-intervention
#' control.
#'
#' @return a list of `twin_intervention` objects.
#' @seealso [intervention()], [scoreInterventions()]
#' @export
#' @examples
#' names(defaultInterventions())
defaultInterventions <- function() {
  list(strengthening   = intervention("strengthening",  scale = list(strength = 1.4),
                                       burden = 0.4, risk = 0.10),
       tone_reduction  = intervention("tone_reduction", scale = list(damping = 0.6, stiffness = 0.7),
                                       burden = 0.3, risk = 0.30),
       stretching      = intervention("stretching",     scale = list(stiffness = 0.6),
                                       burden = 0.2, risk = 0.10),
       orthosis        = intervention("orthosis",       scale = list(stiffness = 1.5),
                                       burden = 0.5, risk = 0.20),
       none            = intervention("none",           burden = 0.0, risk = 0.00))
}

#' Certify a personalised twin before it is allowed to advise (validation gate)
#'
#' The trust gate of the clinical loop. A personalised twin may only inform
#' decisions once it is shown to (1) FIT the patient's data and (2) have
#' IDENTIFIABLE parameters. This recomputes the personalisation residual (fit
#' quality) and a Fisher-information parameter covariance (identifiability), and
#' returns a verdict (`ready` / `provisional` / `not_ready`) with a `[0, 1]`
#' confidence that downstream scoring uses to scale and flag its recommendations.
#'
#' @param fit A `personalized_twin` from [personalizeTwin()].
#' @param imu The IMU data the twin was fitted to.
#' @param dt Sample period (s).
#' @param gyro_var,accel_var Measurement-noise variances (match the fit).
#' @param cv_warn Log-space parameter SD above which a parameter is deemed poorly
#'   identified.
#' @param chi2_ready,chi2_bad Reduced-chi-square thresholds: below `chi2_ready` the
#'   fit is good (a candidate for `ready`); above `chi2_bad` it is a misfit
#'   (`not_ready`).
#' @return a `twin_readiness`: `verdict`, `reduced_chi2` (goodness of fit across
#'   ALL channels), `gyro_nrmse`, the (misfit-inflated) parameter `cov`, per-
#'   parameter `cv`, `confidence`, and the `estimates`/`theta0`/`gyro_bias`.
#' @references Sargent RG (2013) J Simul 7:12-24 (model validation); Bevington &
#'   Robinson (reduced chi-square goodness of fit); Cramer-Rao identifiability.
#' @seealso [personalizeTwin()], [scoreInterventions()], [validateTwin()]
#' @export
#' @examples
#' \donttest{
#' truth <- limbTwin(strength = 7, damping = 0.5)
#' imu <- imuMeasure(simulateTwin(truth, 6, 0.01), imuSensor(gyro_noise = 0.03), seed = 1)
#' fit <- personalizeTwin(imu, limbTwin(strength = 4), estimate = "strength", dt = 0.01)
#' twinReadiness(fit, imu, dt = 0.01)$verdict
#' }
twinReadiness <- function(fit, imu, dt = 0.005, gyro_var = 0.03^2, accel_var = 0.15^2,
                          cv_warn = 0.25, chi2_ready = 3, chi2_bad = 10) {
  stopifnot(inherits(fit, "personalized_twin"))
  est <- fit$estimates; estnm <- names(est); np <- length(est); templ <- fit$twin
  Yg <- imu$gyro; Ya <- imu$ax; Yb <- imu$ay; n <- length(Yg); dur <- (n - 1) * dt
  keep <- floor(0.2 * n):n
  gbias <- fit$gyro_bias %||% 0; theta0 <- fit$theta0 %||% 0
  mk <- function(logp) { tw <- templ; tw[estnm] <- as.list(exp(logp)); tw }
  obj <- function(logp) {           # FULL weighted objective across ALL channels (gyro + accel)
    sim <- simulateTwin(mk(logp), dur, dt, theta0 = theta0, omega0 = Yg[1] - gbias)
    pred <- imuMeasure(sim, noise = FALSE)
    sum((pred$gyro[keep] + gbias - Yg[keep])^2) / gyro_var +
      sum((pred$ax[keep] - Ya[keep])^2 + (pred$ay[keep] - Yb[keep])^2) / accel_var
  }
  lp0 <- log(est)
  # reduced chi-square is the honest goodness-of-fit: ~1 = fits within noise, >>1 =
  # misfit (a wrong local optimum that reproduces one channel is caught by the others)
  dof <- max(1, length(keep) * 3 - (np + 2))
  rchi2 <- obj(lp0) / dof
  sim0 <- simulateTwin(mk(lp0), dur, dt, theta0 = theta0, omega0 = Yg[1] - gbias)
  gyro_nrmse <- sqrt(mean((imuMeasure(sim0, noise = FALSE)$gyro[keep] + gbias - Yg[keep])^2)) /
                  (stats::sd(Yg) + 1e-9)
  fisher <- 0.5 * .num_hessian(obj, lp0)                 # observed Fisher information (log-space)
  cov_raw <- tryCatch(solve(fisher), error = function(e) matrix(Inf, np, np))
  cov <- max(1, rchi2) * (cov_raw + t(cov_raw)) / 2      # inflate by misfit (standard NLS scaling)
  cv <- sqrt(pmax(0, diag(cov))); names(cv) <- estnm     # log-space SD ~ coefficient of variation
  identifiable <- all(is.finite(cv)) && all(cv < cv_warn)
  verdict <- if (!all(is.finite(cv)) || any(diag(cov) < 0) || rchi2 > chi2_bad) "not_ready" else
    if (identifiable && rchi2 < chi2_ready) "ready" else "provisional"
  confidence <- if (!all(is.finite(cv))) 0 else
    max(0, min(1, exp(-max(0, rchi2 - 1) / 2) * exp(-mean(cv))))
  structure(list(verdict = verdict, reduced_chi2 = rchi2, gyro_nrmse = gyro_nrmse,
                 cov = cov, cv = cv, estimates = est, estimated = estnm,
                 confidence = confidence, theta0 = theta0, gyro_bias = gbias),
            class = "twin_readiness")
}

#' @export
print.twin_readiness <- function(x, ...) {
  cat(sprintf("Twin readiness: %s (confidence %.2f, reduced chi-square %.2f)\n",
              toupper(x$verdict), x$confidence, x$reduced_chi2))
  cat(sprintf("  parameter identifiability (log-space SD): %s\n",
              paste(sprintf("%s %.3f", x$estimated, x$cv), collapse = ", ")))
  if (x$verdict != "ready")
    cat("  -> not validated for advice: fit or identifiability inadequate; collect richer data.\n")
  invisible(x)
}

#' Reconcile a mechanistic prediction with external evidence
#'
#' Combines the twin's mechanistic effect estimate with an external (e.g.
#' literature) estimate by inverse-variance (precision) weighting -- the Bayesian
#' conjugate-normal update -- and flags whether the two AGREE. Evidence complements
#' the mechanistic prediction; it is never overridden.
#'
#' @param mech_est,mech_se Mechanistic effect estimate and its SE.
#' @param ev_est,ev_se Evidence-based effect estimate and its SE.
#' @return an `evidence_reconciliation`: combined `estimate`/`se`, the discrepancy
#'   `z`, and `agreement` (`agree`/`disagree`, at |z| < 2).
#' @seealso [scoreInterventions()], [abcCalibration()]
#' @export
#' @examples
#' reconcileEvidence(mech_est = 0.30, mech_se = 0.05, ev_est = 0.26, ev_se = 0.04)
reconcileEvidence <- function(mech_est, mech_se, ev_est, ev_se) {
  wm <- 1 / mech_se^2; we <- 1 / ev_se^2
  est <- (wm * mech_est + we * ev_est) / (wm + we); se <- sqrt(1 / (wm + we))
  z <- abs(mech_est - ev_est) / sqrt(mech_se^2 + ev_se^2)
  structure(list(estimate = est, se = se, z = z,
                 agreement = if (z < 2) "agree" else "disagree",
                 mechanistic = mech_est, evidence = ev_est),
            class = "evidence_reconciliation")
}

#' @export
print.evidence_reconciliation <- function(x, ...) {
  cat(sprintf("Evidence reconciliation: mechanistic %.3f vs evidence %.3f -> %.3f (%s, z=%.1f)\n",
              x$mechanistic, x$evidence, x$estimate, x$agreement, x$z)); invisible(x)
}

#' Score and rank candidate interventions with predictive uncertainty
#'
#' Runs each candidate intervention over the twin's PARAMETER POSTERIOR (sampled
#' from the [twinReadiness()] Fisher covariance), so every candidate gets a
#' predicted effect WITH a predictive interval and the probability that the effect
#' exceeds a minimal clinically important difference (MCID) -- not a point estimate.
#' Candidates are ranked by a multi-criteria score combining efficacy, certainty
#' (probability of benefit), risk and burden, scaled by the twin's readiness
#' confidence; optional per-candidate evidence is reconciled in, and low-confidence
#' or evidence-disagreeing recommendations are flagged.
#'
#' @param twin A personalised `movement_twin` (e.g. `fit$twin`).
#' @param candidates A list of `twin_intervention` (see [defaultInterventions()]).
#' @param readiness A `twin_readiness`; supplies the parameter posterior and the
#'   confidence. If `NULL`, a point estimate is used and the result is marked
#'   `unvalidated`.
#' @param outcome Outcome name (`"rom"`, `"peak_velocity"`, ...).
#' @param mcid Minimal clinically important difference; defaults to 10% of the
#'   baseline outcome.
#' @param n_mc Posterior samples for uncertainty propagation.
#' @param weights Named weights `c(efficacy, certainty, risk, burden)`.
#' @param duration,dt,theta0 Simulation settings; `theta0` defaults to the fit's.
#' @param seed Optional RNG seed.
#' @return an `intervention_scores`: a ranked `ranking` data frame (effect,
#'   `ci_lo`/`ci_hi`, `p_benefit`, risk, burden, reconciled evidence, `score`,
#'   `flag`), the `mcid`, and the `readiness_confidence`/`verdict`.
#' @seealso [twinReadiness()], [reconcileEvidence()], [clinicalDecision()]
#' @export
#' @examples
#' \donttest{
#' truth <- limbTwin(strength = 6, damping = 0.6, stiffness = 3)
#' imu <- imuMeasure(simulateTwin(truth, 6, 0.01), imuSensor(gyro_noise = 0.03), seed = 1)
#' fit <- personalizeTwin(imu, limbTwin(strength = 4), estimate = "strength", dt = 0.01)
#' rdy <- twinReadiness(fit, imu, dt = 0.01)
#' scoreInterventions(fit$twin, defaultInterventions(), rdy, n_mc = 50, seed = 1)$ranking
#' }
scoreInterventions <- function(twin, candidates, readiness = NULL, outcome = "rom",
                               mcid = NULL, n_mc = 200,
                               weights = c(efficacy = 1, certainty = 1, risk = 1, burden = 0.5),
                               duration = 6, dt = 0.005, theta0 = NULL, seed = NULL) {
  stopifnot(inherits(twin, "movement_twin"))
  if (!is.null(seed)) set.seed(seed)
  th0 <- theta0 %||% (if (!is.null(readiness)) readiness$theta0 else 0.1)
  pd <- .posterior_draw(readiness, n_mc)                  # parameter posterior (or point estimate)
  base_out <- .twin_outcomes(simulateTwin(twin, duration, dt, th0, 0))[[outcome]]
  if (is.null(mcid)) mcid <- 0.1 * abs(base_out)
  rows <- lapply(candidates, function(cand) {
    d <- .delta_samples(twin, cand, pd$S, pd$estnm, outcome, duration, dt, th0)
    nmc <- length(d); se <- if (nmc > 1) stats::sd(d) else NA_real_
    rec <- if (!is.null(cand$evidence))
      reconcileEvidence(mean(d), max(se, 1e-6, na.rm = TRUE), cand$evidence$effect, cand$evidence$se) else NULL
    data.frame(intervention = cand$name, effect = mean(d),
               ci_lo = if (nmc > 1) stats::quantile(d, 0.025, names = FALSE) else mean(d),
               ci_hi = if (nmc > 1) stats::quantile(d, 0.975, names = FALSE) else mean(d),
               p_benefit = if (nmc > 1) mean(d > mcid) else as.numeric(mean(d) > mcid),
               risk = cand$risk, burden = cand$burden,
               reconciled = if (is.null(rec)) NA_real_ else rec$estimate,
               agree = if (is.null(rec)) NA_character_ else rec$agreement,
               stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows); rownames(df) <- NULL
  eff <- df$effect
  effn <- if (diff(range(eff)) > 1e-12) (eff - min(eff)) / diff(range(eff)) else rep(0, length(eff))
  conf <- if (!is.null(readiness)) readiness$confidence else 1
  score <- weights["efficacy"] * effn + weights["certainty"] * df$p_benefit -
           weights["risk"] * df$risk - weights["burden"] * df$burden
  score <- as.numeric(score) * conf
  disagree <- !is.na(df$agree) & df$agree == "disagree"
  score[disagree] <- score[disagree] * 0.5               # penalise evidence disagreement
  df$score <- round(score, 3)
  # flag = do NOT trust this number (low twin confidence or evidence conflict) -- NOT
  # the same as "not recommended" (a confidently ineffective option is not flagged)
  df$flag <- ifelse(conf < 0.5 | disagree, "low-confidence", "")
  df <- df[order(-df$score), ]; rownames(df) <- NULL
  numcols <- c("effect", "ci_lo", "ci_hi", "p_benefit", "reconciled", "score")
  df[numcols] <- lapply(df[numcols], function(v) round(v, 4))
  structure(list(ranking = df, outcome = outcome, mcid = mcid,
                 readiness_confidence = conf,
                 verdict = if (!is.null(readiness)) readiness$verdict else "unvalidated"),
            class = "intervention_scores")
}

#' @export
print.intervention_scores <- function(x, ...) {
  cat(sprintf("Intervention ranking (outcome '%s', MCID %.3f, twin %s, confidence %.2f)\n",
              x$outcome, x$mcid, toupper(x$verdict), x$readiness_confidence))
  for (i in seq_len(nrow(x$ranking))) {
    r <- x$ranking[i, ]
    cat(sprintf("  %d. %-14s score %+.3f | %+.3f [%+.3f, %+.3f], P(benefit) %.0f%%%s\n",
                i, r$intervention, r$score, r$effect, r$ci_lo, r$ci_hi, 100 * r$p_benefit,
                if (nzchar(r$flag)) paste0("  [", r$flag, "]") else ""))
  }
  cat("  Research prototype: predictions are conditional on the twin being validated\n",
      "  for this patient; reconcile with clinical judgement and evidence. Not a medical device.\n", sep = "")
  invisible(x)
}

#' Evaluate a delivered intervention against its prediction (close the loop)
#'
#' After an intervention is delivered and a follow-up outcome is measured, check
#' whether the observed change fell within the twin's predicted interval -- the
#' empirical test that closes the personalise -> predict -> act -> evaluate loop
#' and calibrates trust in the twin.
#'
#' @param twin The personalised `movement_twin`.
#' @param chosen The delivered `twin_intervention`.
#' @param observed_delta The measured outcome change at follow-up.
#' @param readiness A `twin_readiness` (for the parameter posterior).
#' @param outcome_sd Outcome observation/session variability (SD). The check uses a
#'   PREDICTION interval combining the parameter (epistemic) uncertainty with this
#'   (aleatoric) noise, so a realistic follow-up is covered -- a parameter-only
#'   interval is over-narrow for a real measurement.
#' @param outcome,n_mc,duration,dt,theta0,seed As in [scoreInterventions()].
#' @return an `intervention_evaluation`: `predicted` effect and prediction `ci`, the
#'   `observed` change, and whether it was `within_interval`.
#' @seealso [scoreInterventions()], [clinicalDecision()]
#' @export
#' @examples
#' \donttest{
#' truth <- limbTwin(strength = 6)
#' imu <- imuMeasure(simulateTwin(truth, 6, 0.01), imuSensor(gyro_noise = 0.03), seed = 1)
#' fit <- personalizeTwin(imu, limbTwin(strength = 4), estimate = "strength", dt = 0.01)
#' rdy <- twinReadiness(fit, imu, dt = 0.01)
#' evaluateIntervention(fit$twin, defaultInterventions()$strengthening, 1.5, rdy,
#'                      outcome_sd = 0.1, n_mc = 50)
#' }
evaluateIntervention <- function(twin, chosen, observed_delta, readiness = NULL,
                                 outcome = "rom", outcome_sd = 0, n_mc = 200,
                                 duration = 6, dt = 0.005, theta0 = NULL, seed = NULL) {
  stopifnot(inherits(twin, "movement_twin"))
  if (!is.null(seed)) set.seed(seed)
  th0 <- theta0 %||% (if (!is.null(readiness)) readiness$theta0 else 0.1)
  pd <- .posterior_draw(readiness, n_mc)
  d <- .delta_samples(twin, chosen, pd$S, pd$estnm, outcome, duration, dt, th0)
  pred <- if (outcome_sd > 0) d + stats::rnorm(length(d), 0, outcome_sd) else d
  ci <- if (length(pred) > 1) stats::quantile(pred, c(0.025, 0.975), names = FALSE)
        else mean(pred) + c(-1.96, 1.96) * outcome_sd
  structure(list(intervention = chosen$name, predicted = mean(d), ci = ci,
                 observed = observed_delta,
                 within_interval = observed_delta >= ci[1] && observed_delta <= ci[2],
                 outcome = outcome), class = "intervention_evaluation")
}

#' @export
print.intervention_evaluation <- function(x, ...) {
  cat(sprintf("Evaluation of '%s' (%s): predicted %+.3f [%+.3f, %+.3f], observed %+.3f -> %s\n",
              x$intervention, x$outcome, x$predicted, x$ci[1], x$ci[2], x$observed,
              if (x$within_interval) "within predicted interval" else "OUTSIDE interval (recalibrate)"))
  invisible(x)
}

#' Run the clinical decision loop end-to-end
#'
#' Personalises a twin to a patient's IMU data, gates it with [twinReadiness()],
#' and -- only reporting through the gate's confidence -- scores and ranks the
#' candidate interventions. A one-call demonstration of the loop; the returned
#' object holds the fit, the readiness verdict and the ranked recommendations.
#'
#' @param imu Patient IMU data (as from [imuMeasure()]).
#' @param prior A `movement_twin` prior/initial guess for [personalizeTwin()].
#' @param estimate Parameters to personalise.
#' @param candidates Candidate interventions (default [defaultInterventions()]).
#' @param outcome Outcome to optimise.
#' @param method Personalisation method (`"optim"`/`"ukf"`).
#' @param dt,n_mc,weights,seed As in [scoreInterventions()].
#' @return a `clinical_decision`: `fit`, `readiness`, and `scores`.
#' @seealso [personalizeTwin()], [twinReadiness()], [scoreInterventions()]
#' @export
#' @examples
#' \donttest{
#' truth <- limbTwin(strength = 6, damping = 0.6, stiffness = 3)
#' imu <- imuMeasure(simulateTwin(truth, 6, 0.01), imuSensor(gyro_noise = 0.03), seed = 1)
#' dec <- clinicalDecision(imu, limbTwin(strength = 4), estimate = "strength",
#'                         dt = 0.01, n_mc = 40, seed = 1)
#' dec$scores$ranking$intervention[1]
#' }
clinicalDecision <- function(imu, prior, estimate = "strength",
                             candidates = defaultInterventions(), outcome = "rom",
                             method = c("optim", "ukf"), dt = 0.005, n_mc = 200,
                             weights = c(efficacy = 1, certainty = 1, risk = 1, burden = 0.5),
                             seed = NULL) {
  fit <- personalizeTwin(imu, prior, estimate = estimate, method = match.arg(method), dt = dt)
  rdy <- twinReadiness(fit, imu, dt = dt)
  scores <- scoreInterventions(fit$twin, candidates, readiness = rdy, outcome = outcome,
                               n_mc = n_mc, weights = weights, dt = dt, seed = seed)
  structure(list(fit = fit, readiness = rdy, scores = scores), class = "clinical_decision")
}

#' @export
print.clinical_decision <- function(x, ...) {
  cat("== Clinical decision loop ==\n")
  print(x$readiness); print(x$scores); invisible(x)
}
