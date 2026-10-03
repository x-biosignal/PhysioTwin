# Validity assessment for a simulation model / digital twin.
#
# Simulation validity is NOT the same as matching one dataset: it splits into
# VERIFICATION (are the equations solved correctly?) and VALIDATION (are they the
# right equations for reality?), and for a twin it adds ASSIMILATION CREDIBILITY
# (does the personalised twin predict data it was not fitted on?). This harness
# runs the applicable checks and returns a structured report.
#
#   L0 verification   -- energy conservation (dissipation-free) + RK4 convergence order
#   L1 identifiability -- can the parameters be recovered in principle? (Fisher / CRB)
#   L2 recovery        -- Monte-Carlo parameter recovery: bias, SD, coverage
#   L3 sensitivity     -- Morris elementary-effects screening (which params matter)
#   L4 prediction      -- cross-condition predictive error + face-validity behaviour
#   L5 domain          -- operating envelope where the model/estimator is valid

# concatenated noise-free IMU output vector for a twin (the "measurement function")
.val_output <- function(tw, duration, dt, theta0) {
  sim <- simulateTwin(tw, duration, dt, theta0 = theta0)
  m <- imuMeasure(sim, noise = FALSE)
  c(m$gyro, m$ax, m$ay)
}

# --- L0: energy conservation (no active drive, no damping = conservative) -----
.val_energy <- function(twin, dt) {
  tw <- twin; tw$strength <- 0; tw$damping <- 0
  sim <- simulateTwin(tw, duration = 6, dt = dt, theta0 = 0.4)
  E <- 0.5 * tw$inertia * sim$omega^2 +
    0.5 * tw$stiffness * (sim$theta - tw$rest)^2 +
    tw$mass * tw$g * tw$com * sin(sim$theta)
  diff(range(E)) / max(mean(abs(E)), .Machine$double.eps)
}

# --- L0: observed order of accuracy of the RK4 integrator ---------------------
.val_order <- function(twin) {
  ref <- simulateTwin(twin, duration = 2, dt = 2.5e-4, theta0 = 0.3)
  end <- function(dt) { s <- simulateTwin(twin, 2, dt, theta0 = 0.3); s$theta[length(s$theta)] }
   e1 <- abs(end(0.02) - ref$theta[length(ref$theta)])
  e2 <- abs(end(0.01) - ref$theta[length(ref$theta)])
  if (e2 <= 0) return(Inf)
  log2(e1 / e2)                                        # ~4 for RK4
}

# --- L1: local identifiability via the sensitivity / Fisher matrix ------------
.val_identifiability <- function(twin, estimate, duration, dt, theta0, sigma) {
  y0 <- .val_output(twin, duration, dt, theta0)
  S <- vapply(estimate, function(p) {                 # d output / d log(param)
    tw <- twin; h <- 0.01
    tw[[p]] <- twin[[p]] * exp(h); yp <- .val_output(tw, duration, dt, theta0)
    tw[[p]] <- twin[[p]] * exp(-h); ym <- .val_output(tw, duration, dt, theta0)
    (yp - ym) / (2 * h)
  }, numeric(length(y0)))
  FIM <- crossprod(S) / sigma^2
  sv <- svd(S)$d
  list(rank = sum(sv > 1e-8 * max(sv)), n_params = length(estimate),
       condition = if (min(sv) > 0) max(sv) / min(sv) else Inf,
       cramer_rao_sd = sqrt(diag(solve(FIM + diag(1e-12, ncol(FIM))))))
}

# --- L2: Monte-Carlo parameter recovery --------------------------------------
.val_recovery <- function(twin, estimate, duration, dt, theta0, sensor, n_mc, method) {
  truth <- unlist(twin[estimate])
  prior <- twin; for (p in estimate) prior[[p]] <- twin[[p]] * 0.6   # off-by-40% prior
  est <- matrix(NA_real_, n_mc, length(estimate), dimnames = list(NULL, estimate))
  sim <- simulateTwin(twin, duration, dt, theta0 = theta0)
  for (i in seq_len(n_mc)) {
    imu <- imuMeasure(sim, sensor, seed = 1000 + i)
    fit <- personalizeTwin(imu, prior, estimate = estimate, method = method, dt = dt,
                           gyro_var = sensor$gyro_noise^2, accel_var = sensor$accel_noise^2,
                           n_starts = 1L)                 # off-by-40% prior is near truth; keep the MC loop cheap
    est[i, ] <- fit$estimates[estimate]
  }
  bias <- colMeans(est) - truth; sdv <- apply(est, 2, stats::sd)
  cover <- vapply(estimate, function(p)
    mean(abs(est[, p] - truth[p]) <= 0.15 * truth[p]), numeric(1))   # within 15%
  data.frame(param = estimate, true = truth, bias = bias, sd = sdv,
             coverage_15pct = cover, row.names = NULL)
}

# --- L3: Morris elementary-effects sensitivity of an output metric ------------
.val_sensitivity <- function(twin, estimate, duration, dt, theta0, n_traj) {
  metric <- function(tw) diff(range(simulateTwin(tw, duration, dt, theta0 = theta0)$theta))
  ee <- matrix(NA_real_, n_traj, length(estimate), dimnames = list(NULL, estimate))
  for (r in seq_len(n_traj)) {
    base <- twin
    for (p in estimate) base[[p]] <- twin[[p]] * stats::runif(1, 0.7, 1.3)
    m0 <- metric(base)
    for (p in estimate) {
      tw <- base; delta <- 0.1 * twin[[p]]; tw[[p]] <- base[[p]] + delta
      ee[r, p] <- (metric(tw) - m0) / delta
    }
  }
  colMeans(abs(ee))                                    # mu* (mean absolute elementary effect)
}

# --- L4: cross-condition prediction + face validity --------------------------
.val_prediction <- function(twin, estimate, dt, theta0, sensor) {
  fit_tw <- twin; fit_tw$freq <- twin$freq              # condition A = the model freq
  simA <- simulateTwin(fit_tw, 8, dt, theta0 = theta0)
  imuA <- imuMeasure(simA, sensor, seed = 7)
  prior <- twin; for (p in estimate) prior[[p]] <- twin[[p]] * 0.6
  fit <- personalizeTwin(imuA, prior, estimate = estimate, method = "optim", dt = dt,
                         gyro_var = sensor$gyro_noise^2, accel_var = sensor$accel_noise^2,
                         n_starts = 1L)
  condB <- twin; condB$freq <- twin$freq * 1.5          # UNSEEN condition B
  predB <- fit$twin; predB$freq <- twin$freq * 1.5
  yB_true <- .val_output(condB, 8, dt, theta0); yB_pred <- .val_output(predB, 8, dt, theta0)
  rmse <- sqrt(mean((yB_true - yB_pred)^2)) / sqrt(mean(yB_true^2))   # relative
  # face validity: known monotone behaviours
  rom <- function(tw) diff(range(simulateTwin(tw, 8, dt, theta0 = theta0)$theta))
  pv  <- function(tw) max(abs(simulateTwin(tw, 8, dt, theta0 = theta0)$omega))
  up <- function(p, f) { a <- twin; a[[p]] <- twin[[p]] * 1.4; f(a) }
  face <- c(strength_raises_rom = up("strength", rom) > rom(twin),
            damping_lowers_peakvel = up("damping", pv) < pv(twin))
  list(cross_condition_rel_rmse = rmse, face_validity = face)
}

#' Validity assessment of a movement twin (verification & validation)
#'
#' Runs the applicable simulation-validity checks and returns a structured
#' report: verification (energy conservation, integrator order), identifiability
#' (can parameters be recovered in principle), Monte-Carlo recovery (bias / SD /
#' coverage), sensitivity screening (which parameters matter), cross-condition
#' prediction + face validity, and the domain of applicability. This is
#' verification and internal validation on the model itself; validation against a
#' subject's REAL data is done separately by personalising and predicting held-out
#' recordings.
#'
#' @param twin A `movement_twin` (the ground-truth model to assess).
#' @param estimate Parameters whose identifiability/recovery are assessed.
#' @param sensor The `imu_sensor` used for the recovery / prediction checks.
#' @param duration,dt,theta0 Simulation settings for the checks.
#' @param n_mc,n_traj Monte-Carlo recovery repetitions and Morris trajectories.
#' @param method Assimilation method for the recovery check (`"optim"`/`"ukf"`).
#' @param checks Which layers to run (subset of `c("verification",
#'   "identifiability", "recovery", "sensitivity", "prediction", "domain")`).
#' @return a `twin_validity` report.
#' @references Sargent RG (2013) J Simulation 7:12-24; ASME V&V 10/40;
#'   Oberkampf & Roy (2010) Verification and Validation in Scientific Computing.
#' @seealso [personalizeTwin()], [insilicoIntervention()]
#' @export
#' @examples
#' \donttest{
#' rep <- validateTwin(limbTwin(strength = 6, damping = 0.5, stiffness = 6),
#'                     estimate = c("strength", "damping"), n_mc = 8, n_traj = 5)
#' rep$verification
#' }
validateTwin <- function(twin, estimate = c("strength", "damping"),
                         sensor = imuSensor(gyro_noise = 0.03, accel_noise = 0.12),
                         duration = 6, dt = 0.02, theta0 = 0.2, n_mc = 15L, n_traj = 8L,
                         method = "optim",
                         checks = c("verification", "identifiability", "recovery",
                                    "sensitivity", "prediction", "domain")) {
  stopifnot(inherits(twin, "movement_twin"))
  out <- list(twin = twin, estimate = estimate)
  if ("verification" %in% checks)
    out$verification <- list(energy_drift = .val_energy(twin, dt),
                             rk4_order = .val_order(twin))
  if ("identifiability" %in% checks)
    out$identifiability <- .val_identifiability(twin, estimate, duration, dt, theta0,
                                                sensor$gyro_noise)
  if ("recovery" %in% checks)
    out$recovery <- .val_recovery(twin, estimate, duration, dt, theta0, sensor, n_mc, method)
  if ("sensitivity" %in% checks)
    out$sensitivity <- .val_sensitivity(twin, estimate, duration, dt, theta0, n_traj)
  if ("prediction" %in% checks)
    out$prediction <- .val_prediction(twin, estimate, dt, theta0, sensor)
  if ("domain" %in% checks) {
    amp <- diff(range(simulateTwin(twin, duration, dt, theta0 = theta0)$theta))
    out$domain <- list(amplitude_rad = amp,
                       regime = if (amp > 1) "large-amplitude (strongly nonlinear)" else "moderate",
                       recommended_method = if (amp > 1) "optim" else "optim or ukf")
  }
  structure(out, class = "twin_validity")
}

#' @export
print.twin_validity <- function(x, ...) {
  cat("Twin validity report\n")
  if (!is.null(x$verification))
    cat(sprintf("  L0 verification : energy drift %.2e, RK4 order %.2f\n",
                x$verification$energy_drift, x$verification$rk4_order))
  if (!is.null(x$identifiability))
    cat(sprintf("  L1 identifiab.  : rank %d/%d, condition %.1f\n",
                x$identifiability$rank, x$identifiability$n_params, x$identifiability$condition))
  if (!is.null(x$recovery)) {
    r <- x$recovery
    cat("  L2 recovery     :",
        paste(sprintf("%s bias %+.2f cover %.0f%%", r$param, r$bias, 100 * r$coverage_15pct),
              collapse = "; "), "\n")
  }
  if (!is.null(x$sensitivity))
    cat("  L3 sensitivity  :",
        paste(sprintf("%s=%.2f", names(x$sensitivity), x$sensitivity), collapse = " "), "\n")
  if (!is.null(x$prediction))
    cat(sprintf("  L4 prediction   : cross-condition rel-RMSE %.3f, face-validity %d/%d\n",
                x$prediction$cross_condition_rel_rmse, sum(x$prediction$face_validity),
                length(x$prediction$face_validity)))
  if (!is.null(x$domain))
    cat(sprintf("  L5 domain       : amplitude %.2f rad (%s)\n",
                x$domain$amplitude_rad, x$domain$regime))
  invisible(x)
}
