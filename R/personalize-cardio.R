# Layer 3/4 -- multi-system personalisation: fit the closed-loop cardiovascular-
# respiratory twin to a subject's heart-rate and blood-pressure variability.
#
# The open-loop movement twin is personalised from an IMU stream; the closed-loop
# CV-respiratory twin is personalised from the variability SUMMARIES a clinician
# actually measures -- the baroreflex sensitivity, the Mayer-wave (low-frequency)
# power of blood pressure, and the respiratory (high-frequency) power of heart
# rate. Because the map from parameters to these summaries is stochastic and has
# no tractable likelihood, the fit is done by GP-accelerated calibration
# (Bayesian optimisation): a handful of closed-loop simulations locate the
# baroreflex gains and respiratory modulation that reproduce the subject's HRV/BPV.
# Dependency-free base R; uses gpCalibrate() and cardioRespiratory().

# the HRV/BPV summaries a subject is characterised by
.cr_summary <- function(cr) c(brs = cr$brs, bp_lf = cr$bpv$lf_power, hr_hf = cr$hrv$hf_power)

#' Summarise a closed-loop cardiovascular-respiratory recording
#'
#' The heart-rate/blood-pressure variability summaries used to personalise the
#' closed-loop twin: baroreflex sensitivity (`brs`, ms/mmHg), the Mayer-wave
#' low-frequency power of blood pressure (`bp_lf`) and the respiratory
#' high-frequency power of heart rate (`hr_hf`).
#'
#' @param cr a `cardio_respiratory` object from [cardioRespiratory()].
#' @return a named numeric vector `c(brs, bp_lf, hr_hf)`.
#' @seealso [cardioRespiratory()], [personalizeCardio()]
#' @export
#' @examples
#' cardioSummary(cardioRespiratory(duration = 120, seed = 1))
cardioSummary <- function(cr) {
  stopifnot(inherits(cr, "cardio_respiratory")); .cr_summary(cr)
}

#' Personalise the closed-loop CV-respiratory twin to real HRV/BPV summaries
#'
#' Fits the closed-loop [cardioRespiratory()] model's baroreflex gains and
#' respiratory modulation to a subject's measured variability summaries (baroreflex
#' sensitivity, blood-pressure Mayer power, heart-rate respiratory power) by
#' GP-accelerated calibration ([gpCalibrate()]). Because the summaries have very
#' different scales, the discrepancy is relative. The returned twin, driven with
#' the fitted parameters, reproduces the subject's HRV/BPV.
#'
#' @param target Named numeric vector of the observed summaries to match; any of
#'   `brs`, `bp_lf`, `hr_hf` (see [cardioSummary()]).
#' @param ranges Named list of `c(lower, upper)` bounds for the parameters to fit
#'   (default `g_vagal`, `g_resist`, `rsa_amp`).
#' @param fixed Named list of other [cardioRespiratory()] parameters held fixed.
#' @param duration Simulation length per evaluation (s).
#' @param n_init,n_iter Bayesian-optimisation initial design and iterations.
#' @param n_rep Simulations averaged per evaluation to stabilise the (stochastic)
#'   variability summaries -- the blood-pressure Mayer power is noisy, so averaging
#'   makes the inverse problem well-posed.
#' @param seed RNG seed (fixes the stochastic loss so the calibration is deterministic).
#' @return a `cardio_personalization`: the fitted `estimates`, the personalised
#'   `twin`, the `target` and `achieved` summaries, the relative error `rel_error`,
#'   and the `calibration` object.
#' @seealso [cardioRespiratory()], [cardioSummary()], [gpCalibrate()], [personalizeTwin()]
#' @export
#' @examples
#' \donttest{
#' subject <- cardioSummary(cardioRespiratory(duration = 240, g_vagal = 0.013,
#'                          g_resist = 0.008, rsa_amp = 0.045, seed = 3))
#' fit <- personalizeCardio(subject, n_init = 8, n_iter = 15, seed = 1)
#' fit$achieved            # reproduces the subject's brs / bp_lf / hr_hf
#' }
personalizeCardio <- function(target,
                              ranges = list(g_vagal = c(0.004, 0.020),
                                            g_resist = c(0.003, 0.020),
                                            rsa_amp = c(0.010, 0.060)),
                              fixed = list(), duration = 300, n_init = 12, n_iter = 25,
                              n_rep = 3, seed = 1) {
  tnm <- names(target)
  if (is.null(tnm) || !all(tnm %in% c("brs", "bp_lf", "hr_hf")))
    stop("target must be a named vector of summaries (any of brs, bp_lf, hr_hf).", call. = FALSE)
  scale <- pmax(abs(target), 1e-6); seeds <- seed + seq_len(n_rep)
  sim_summary <- function(p) {          # mean summaries over n_rep seeds (stabilise the noisy powers)
    old <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
    S <- vapply(seeds, function(sd) {
      cr <- do.call(cardioRespiratory, c(as.list(p), fixed, list(duration = duration, seed = sd)))
      .cr_summary(cr)[tnm]
    }, numeric(length(tnm)))
    if (!is.null(old)) assign(".Random.seed", old, envir = .GlobalEnv)  # seeded sims must not clobber the BO stream
    rowMeans(matrix(S, nrow = length(tnm), dimnames = list(tnm, NULL)))
  }
  loss <- function(p) sum(((sim_summary(p) - target) / scale)^2)
  cal <- gpCalibrate(loss, ranges, n_init = n_init, n_iter = n_iter, seed = seed)
  lo <- vapply(ranges, `[`, numeric(1), 1); hi <- vapply(ranges, `[`, numeric(1), 2)
  op <- tryCatch(stats::optim(cal$par, loss, method = "L-BFGS-B", lower = lo, upper = hi,
                              control = list(maxit = 25)), error = function(e) NULL)  # local polish
  best <- if (!is.null(op) && op$value < cal$value) op$par else cal$par
  cal$par <- best
  ach <- sim_summary(best)
  cr <- do.call(cardioRespiratory, c(as.list(best), fixed, list(duration = duration, seed = seeds[1])))
  structure(list(estimates = cal$par, twin = cr, target = target, achieved = ach,
                 rel_error = (ach - target) / scale, calibration = cal),
            class = "cardio_personalization")
}

#' @export
print.cardio_personalization <- function(x, ...) {
  cat("Cardio-respiratory personalisation:\n")
  cat(sprintf("  estimates: %s\n",
              paste(sprintf("%s=%.4g", names(x$estimates), x$estimates), collapse = ", ")))
  for (nm in names(x$target))
    cat(sprintf("  %-6s target %.3g -> achieved %.3g (%+.0f%%)\n",
                nm, x$target[[nm]], x$achieved[[nm]], 100 * x$rel_error[[nm]]))
  invisible(x)
}
