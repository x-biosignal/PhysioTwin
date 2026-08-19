# Bayesian optimal experimental design for identifying twin parameters.
#
# Before collecting data, WHICH experimental condition should be run? The one that
# will best IDENTIFY the parameters -- that is, the condition whose data carry the
# most Fisher information about them, so the estimator's variance (the Cramer-Rao
# lower bound) is smallest. optimalDesign() scores a grid of candidate conditions
# (e.g. driving frequencies) by that information and returns the most informative,
# turning the twin into a design tool: run the experiment worth running. Mirrors
# the Fisher/sensitivity construction of the validity harness (validate.R).

# Fisher information of the (inverse-noise-weighted) IMU observable about the
# parameters `estimate`, at a given twin/design. Channels are weighted by their own
# noise, so the gyro and accelerometer contribute in proportion to their precision.
.design_fim <- function(twin, estimate, duration, dt, theta0, gyro_sd, accel_sd) {
  y0 <- .val_output(twin, duration, dt, theta0); n <- length(y0) / 3
  w <- c(rep(1 / gyro_sd, n), rep(1 / accel_sd, 2 * n))     # inverse-noise weights per channel
  S <- vapply(estimate, function(p) {                        # d output / d log(param), central diff
    tw <- twin; h <- 0.01
    tw[[p]] <- twin[[p]] * exp(h);  yp <- .val_output(tw, duration, dt, theta0)
    tw[[p]] <- twin[[p]] * exp(-h); ym <- .val_output(tw, duration, dt, theta0)
    (yp - ym) / (2 * h)
  }, numeric(length(y0)))
  crossprod(S * w)                                           # FIM = t(Sw) %*% Sw
}

#' Bayesian optimal experimental design for twin identification
#'
#' Chooses the experimental condition that will best identify a movement twin's
#' parameters. For each candidate value of a design variable (e.g. the driving
#' frequency), it forms the Fisher information matrix of the noise-weighted IMU
#' observable about the parameters and scores the design by an optimality criterion
#' -- D-optimality (the log-determinant of the information, maximising overall
#' information) or A-optimality (minimising the total Cramer-Rao variance). The most
#' informative condition is the experiment worth running to pin the parameters down.
#'
#' @param twin A `movement_twin` giving the assumed (prior) parameter values.
#' @param estimate Parameter(s) whose identifiability is optimised.
#' @param designs Numeric vector of candidate values for the design variable.
#' @param design_var The twin field the design varies (e.g. `"freq"`, `"duration"`).
#' @param sensor An `imuSensor` (its `gyro_noise`/`accel_noise` weight the channels).
#' @param duration,dt,theta0 Simulation settings for the observable.
#' @param criterion `"D"` (max log-det Fisher information) or `"A"` (min total CRLB variance).
#' @return an `optimal_design`: a `table` (`design`, `logdet_fim`, `score`, and a
#'   `crlb_<param>` Cramer-Rao SD per estimated parameter), the `best` design, and
#'   the `estimate`/`design_var`/`criterion`.
#' @references Chaloner K, Verdinelli I (1995) Bayesian experimental design: a
#'   review. Statist Sci 10:273-304; Fedorov (1972) D-optimality.
#' @seealso [personalizeTwin()], [validateTwin()], [profileLikelihood()]
#' @export
#' @examples
#' od <- optimalDesign(limbTwin(), estimate = "strength",
#'                     designs = seq(0.4, 2.4, by = 0.2), design_var = "freq",
#'                     duration = 5, dt = 0.02)
#' od$best        # the most informative driving frequency
optimalDesign <- function(twin, estimate = "strength", designs, design_var = "freq",
                          sensor = imuSensor(), duration = 6, dt = 0.01, theta0 = 0.1,
                          criterion = c("D", "A")) {
  stopifnot(inherits(twin, "movement_twin"))
  criterion <- match.arg(criterion); estimate <- as.character(estimate)
  gsd <- sensor$gyro_noise; asd <- sensor$accel_noise
  rows <- lapply(designs, function(d) {
    tw <- twin; tw[[design_var]] <- d
    FIM <- .design_fim(tw, estimate, duration, dt, theta0, gsd, asd)
    crb <- tryCatch(sqrt(diag(solve(FIM + diag(1e-12, ncol(FIM))))),
                    error = function(e) rep(Inf, length(estimate)))
    ld <- as.numeric(determinant(FIM, logarithm = TRUE)$modulus)
    score <- if (criterion == "D") ld else -sum(crb^2)      # D: log-det FIM; A: -total CRLB variance
    data.frame(design = d, logdet_fim = ld, score = score,
               stats::setNames(as.list(crb), paste0("crlb_", estimate)), check.names = FALSE)
  })
  tab <- do.call(rbind, rows); rownames(tab) <- NULL
  structure(list(table = tab, best = tab$design[which.max(tab$score)],
                 estimate = estimate, design_var = design_var, criterion = criterion),
            class = "optimal_design")
}

#' @export
print.optimal_design <- function(x, ...) {
  cat(sprintf("Optimal experimental design (%s-optimal) -- identify %s\n",
              x$criterion, paste(x$estimate, collapse = ", ")))
  cat(sprintf("  best %s = %.3g of %d candidates (log-det FIM %.2f)\n",
              x$design_var, x$best, nrow(x$table),
              x$table$logdet_fim[which.max(x$table$score)]))
  invisible(x)
}
