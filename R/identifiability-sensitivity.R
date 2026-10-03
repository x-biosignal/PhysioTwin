# Deeper validity tools: practical identifiability (profile likelihood) and
# variance-based global sensitivity (Sobol indices).
#
# The Fisher/Cramer-Rao check in validateTwin() is LOCAL; these are the global
# tools. A profile likelihood fixes one parameter across a grid, re-optimising the
# rest, and traces the resulting fit cost -- a well-defined minimum with a
# threshold-based confidence interval means the parameter is PRACTICALLY
# identifiable; a flat profile means it is not. Sobol indices apportion the
# variance of a model output among the parameters (first-order and total-effect),
# the definitive answer to "which parameters matter".

# IMU-fit cost of a twin against a recording (same objective as personalizeTwin optim)
.id_cost <- function(tw, imu, dt, gyro_var, accel_var, theta0, gbias) {
  n <- nrow(imu); dur <- (n - 1) * dt; keep <- floor(0.2 * n):n
  sim <- simulateTwin(tw, dur, dt, theta0 = theta0, omega0 = imu$gyro[1] - gbias)
  pred <- imuMeasure(sim, noise = FALSE)
  sum(((pred$gyro[keep] + gbias) - imu$gyro[keep])^2) / gyro_var +
    sum((pred$ax[keep] - imu$ax[keep])^2 + (pred$ay[keep] - imu$ay[keep])^2) / accel_var
}

#' Profile-likelihood practical identifiability
#'
#' Assesses whether a parameter is *practically* identifiable from a recording:
#' it is fixed across a grid of values and the remaining `nuisance` parameters
#' (plus a start angle and gyro bias) are re-optimised, tracing the fit cost. A
#' clear minimum bounded by the chi-square threshold gives a confidence interval;
#' a flat profile means the parameter cannot be pinned down.
#'
#' @param imu An IMU recording (`gyro, ax, ay`), real or simulated.
#' @param prior A `movement_twin` giving the fixed parameters and start point.
#' @param param The parameter to profile.
#' @param over Length-2 multiplicative range around `prior[[param]]` (default
#'   `c(0.4, 2.0)`), or an explicit numeric grid.
#' @param nuisance Other parameters re-optimised at each grid point.
#' @param n_grid Grid size (default 15).
#' @param dt,gyro_var,accel_var Sample period and measurement variances.
#' @return a `profile_likelihood`: the `profile` (grid value, cost), the cost
#'   `minimum`, the 95% `ci` (or `NA` if open), and `identifiable` (logical).
#' @references Raue A, et al. (2009) Bioinformatics 25:1923-1929.
#' @seealso [validateTwin()], [sobolIndices()]
#' @export
#' @examples
#' \donttest{
#' tw <- limbTwin(strength = 6, damping = 0.5, stiffness = 6)
#' imu <- imuMeasure(simulateTwin(tw, 8, 0.02, theta0 = 0.2),
#'                   imuSensor(gyro_noise = 0.03, accel_noise = 0.12), seed = 1)
#' profileLikelihood(imu, tw, "strength", nuisance = "damping", n_grid = 9)$ci
#' }
profileLikelihood <- function(imu, prior, param, over = c(0.4, 2.0),
                              nuisance = character(0), n_grid = 15L, dt = 0.02,
                              gyro_var = 0.03^2, accel_var = 0.12^2) {
  stopifnot(inherits(prior, "movement_twin"))
  grid <- if (length(over) == 2) prior[[param]] * seq(over[1], over[2], length.out = n_grid) else over
  cost_at <- function(val) {
    base <- prior; base[[param]] <- val
    if (length(nuisance) == 0) return(
      stats::optim(c(0, 0), function(p) { b <- base
        .id_cost(b, imu, dt, gyro_var, accel_var, p[1], p[2]) }, method = "Nelder-Mead")$value)
    obj <- function(p) { b <- base; b[nuisance] <- as.list(exp(p[seq_along(nuisance)]))
      .id_cost(b, imu, dt, gyro_var, accel_var, p[length(nuisance) + 1], p[length(nuisance) + 2]) }
    stats::optim(c(log(unlist(prior[nuisance])), 0, 0), obj, method = "Nelder-Mead",
                 control = list(maxit = 400))$value
  }
  cost <- vapply(grid, cost_at, numeric(1))
  imin <- which.min(cost); mn <- cost[imin]; thr <- mn + 3.84   # chi-square(1, 0.95)
  interp <- function(i) grid[i] + (thr - cost[i]) / (cost[i + 1] - cost[i]) * (grid[i + 1] - grid[i])
  lo <- NA_real_
  if (imin > 1) for (i in rev(seq_len(imin - 1))) if (cost[i] >= thr) { lo <- interp(i); break }
  hi <- NA_real_
  if (imin < length(grid)) for (i in seq(imin, length(grid) - 1)) if (cost[i + 1] >= thr) { hi <- interp(i); break }
  ci <- c(lo, hi)
  # identifiable = the fit cost rises past the threshold on BOTH sides within the
  # explored range (a bounded minimum); a flat/open profile is not identifiable.
  structure(list(profile = data.frame(value = grid, cost = cost), minimum = mn,
                 best = grid[imin], ci = ci, identifiable = all(is.finite(ci)),
                 param = param), class = "profile_likelihood")
}

#' @export
print.profile_likelihood <- function(x, ...) {
  cat(sprintf("Profile likelihood [%s] -- %s", x$param,
              if (x$identifiable) sprintf("identifiable, 95%% CI [%.2f, %.2f]\n", x$ci[1], x$ci[2])
              else "NOT practically identifiable (open profile)\n"))
  invisible(x)
}

#' Variance-based global sensitivity (Sobol indices)
#'
#' Apportions the variance of a model output among the parameters using the
#' Saltelli/Jansen estimator: first-order indices `Si` (each parameter's own
#' contribution) and total-effect indices `STi` (including interactions). A near-
#' zero total index means the parameter is non-influential (and typically
#' unidentifiable).
#'
#' @param ranges A named list of length-2 parameter ranges (uniform sampling).
#' @param output A function `function(twin)` returning a scalar output metric;
#'   default the movement range of motion.
#' @param base A `movement_twin` for the parameters NOT in `ranges`.
#' @param n Base sample size (total model evaluations `= n * (k + 2)`; default 256).
#' @param duration,dt,theta0 Simulation settings for the default output.
#' @param seed Optional RNG seed.
#' @return a `sobol_indices` data frame: `param`, first-order `Si`, total-effect
#'   `STi`.
#' @references Saltelli A, et al. (2010) Comput Phys Commun 181:259-270;
#'   Jansen MJW (1999).
#' @seealso [validateTwin()], [profileLikelihood()]
#' @export
#' @examples
#' sobolIndices(list(strength = c(3, 9), damping = c(0.2, 1.0), stiffness = c(2, 6)),
#'              n = 64, seed = 1)
sobolIndices <- function(ranges, output = NULL, base = limbTwin(), n = 256L,
                         duration = 6, dt = 0.02, theta0 = 0.2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  k <- length(ranges); nm <- names(ranges)
  if (is.null(output)) output <- function(tw) diff(range(simulateTwin(tw, duration, dt, theta0 = theta0)$theta))
  samp <- function() vapply(nm, function(p) stats::runif(n, ranges[[p]][1], ranges[[p]][2]), numeric(n))
  A <- samp(); B <- samp()
  ev <- function(M) vapply(seq_len(n), function(i) {
    tw <- base; tw[nm] <- as.list(M[i, ]); output(tw) }, numeric(1))
  yA <- ev(A); yB <- ev(B); vY <- stats::var(c(yA, yB))
  Si <- STi <- numeric(k)
  for (j in seq_len(k)) {
    AB <- A; AB[, j] <- B[, j]; yAB <- ev(AB)
    Si[j]  <- mean(yB * (yAB - yA)) / vY                # Saltelli first-order
    STi[j] <- mean((yA - yAB)^2) / (2 * vY)             # Jansen total-effect
  }
  structure(data.frame(param = nm, Si = pmax(pmin(Si, 1), 0), STi = pmax(pmin(STi, 1), 0),
                       row.names = NULL), class = c("sobol_indices", "data.frame"))
}

#' @export
print.sobol_indices <- function(x, ...) {
  cat("Sobol sensitivity indices (Si first-order, STi total-effect)\n")
  for (i in seq_len(nrow(x)))
    cat(sprintf("  %-12s Si=%.2f  STi=%.2f\n", x$param[i], x$Si[i], x$STi[i]))
  invisible(x)
}
