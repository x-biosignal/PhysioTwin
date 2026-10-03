# Layer 3 -- data assimilation: likelihood-free Bayesian calibration (ABC).
#
# Mechanistic twins are easy to SIMULATE but their likelihood is usually
# intractable. Approximate Bayesian Computation sidesteps the likelihood
# entirely: draw parameters from the prior, simulate, keep the draws whose
# simulated summary statistics land closest to the observed ones. The kept draws
# approximate the posterior -- full uncertainty over the mechanistic parameters,
# not just a point estimate. An optional local-linear regression adjustment
# (Beaumont et al. 2002) sharpens the approximation.
# Dependency-free base R.

#' Approximate Bayesian Computation (likelihood-free parameter calibration)
#'
#' Calibrates a simulator's parameters to observed data without a likelihood:
#' parameters are drawn from `prior`, passed through `simulate` to produce summary
#' statistics, and the draws whose summaries fall closest to `observed` (the
#' smallest-distance `accept` fraction) are retained as posterior samples.
#' Distances are Euclidean on summaries scaled by their simulated SD, so
#' statistics on different scales contribute comparably. With `regression = TRUE`
#' the accepted draws are corrected by a distance-weighted local-linear regression
#' onto the observed summaries.
#'
#' @param simulate `function(theta)` returning a numeric vector of summary
#'   statistics for a parameter vector `theta`.
#' @param observed Observed summary statistics (same length as `simulate`'s output).
#' @param prior `function()` returning one prior draw (a named numeric vector).
#' @param n Number of prior simulations.
#' @param accept Fraction of draws to retain (the ABC tolerance as a quantile).
#' @param regression If `TRUE`, apply the Beaumont local-linear regression adjustment.
#' @param seed Optional RNG seed.
#' @return an `abc_posterior`: `samples` (accepted, adjusted parameter draws),
#'   `raw` (accepted draws before adjustment), `mean`, `sd`, `quantiles` (2.5/50/
#'   97.5%), `threshold` (accepted distance) and `accept_rate`.
#' @references Beaumont MA, Zhang W, Balding DJ (2002) Genetics 162:2025-2035;
#'   Pritchard JK et al. (1999).
#' @seealso [personalizeTwin()], [profileLikelihood()], [gpEmulator()]
#' @export
#' @examples
#' # recover the mean/SD of a Gaussian from its sample mean and SD
#' set.seed(1); obs <- c(mean = 2, sd = 1.5)
#' sim <- function(th) { x <- rnorm(200, th[1], th[2]); c(mean(x), sd(x)) }
#' pri <- function() c(m = runif(1, -3, 7), s = runif(1, 0.2, 4))
#' post <- abcCalibration(sim, obs, pri, n = 2000, accept = 0.02, seed = 1)
#' post$mean
abcCalibration <- function(simulate, observed, prior, n = 5000, accept = 0.02,
                           regression = TRUE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  observed <- as.numeric(observed)
  th1 <- prior(); p <- length(th1); q <- length(observed)
  pnames <- names(th1); if (is.null(pnames)) pnames <- paste0("theta", seq_len(p))
  P <- matrix(0, n, p, dimnames = list(NULL, pnames)); S <- matrix(0, n, q)
  P[1, ] <- th1; S[1, ] <- as.numeric(simulate(th1))
  for (i in 2:n) { th <- prior(); P[i, ] <- th; S[i, ] <- as.numeric(simulate(th)) }
  sdev <- apply(S, 2, stats::sd); sdev[sdev < 1e-8] <- 1
  Z <- sweep(sweep(S, 2, observed, "-"), 2, sdev, "/")     # scaled summary deviations
  D <- sqrt(rowSums(Z^2))
  thr <- stats::quantile(D, accept, names = FALSE)
  keep <- which(D <= thr)
  Pa <- P[keep, , drop = FALSE]; Zk <- Z[keep, , drop = FALSE]; Dk <- D[keep]
  Padj <- Pa
  if (regression && length(keep) > (q + 2)) {              # Beaumont local-linear correction
    w <- 1 - (Dk / max(Dk))^2                              # Epanechnikov weights
    for (j in seq_len(p)) {
      b <- stats::lm.wfit(cbind(1, Zk), Pa[, j], w)$coefficients
      b[is.na(b)] <- 0
      Padj[, j] <- Pa[, j] - Zk %*% b[-1]                  # value extrapolated to Z = 0 (observed)
    }
  }
  qs <- apply(Padj, 2, stats::quantile, probs = c(0.025, 0.5, 0.975), names = FALSE)
  structure(list(samples = Padj, raw = Pa,
                 mean = colMeans(Padj), sd = apply(Padj, 2, stats::sd),
                 quantiles = qs, threshold = thr, accept_rate = length(keep) / n,
                 pnames = pnames), class = "abc_posterior")
}

#' @export
print.abc_posterior <- function(x, ...) {
  cat(sprintf("ABC posterior -- %d accepted (%.1f%%):\n", nrow(x$samples), 100 * x$accept_rate))
  for (j in seq_along(x$pnames))
    cat(sprintf("  %-10s %.3f  [%.3f, %.3f]\n", x$pnames[j], x$mean[j],
                x$quantiles[1, j], x$quantiles[3, j]))
  invisible(x)
}
