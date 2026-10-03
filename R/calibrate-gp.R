# Layer 3 -- data assimilation: GP-accelerated calibration (Bayesian optimisation).
#
# When each evaluation of the fit loss is an expensive simulation, calibrating by
# grid or Nelder-Mead wastes runs. Bayesian optimisation fits a Gaussian-process
# surrogate to the (parameters -> loss) points seen so far and uses it to choose
# the NEXT parameters to try -- the point of maximum EXPECTED IMPROVEMENT, which
# balances exploiting the current best against exploring where the surrogate is
# uncertain. It reaches a good calibration in a handful of simulations. Uses
# gpEmulator(); dependency-free base R.

#' GP-accelerated calibration (Bayesian optimisation of an expensive loss)
#'
#' Minimises an expensive loss over a bounded parameter box by Bayesian
#' optimisation: it fits a [gpEmulator()] surrogate to the evaluated
#' (parameters, loss) points and repeatedly evaluates the loss at the point of
#' maximum expected improvement. Far fewer loss evaluations than grid search or a
#' local optimiser -- the tool for calibrating a slow simulator to data.
#'
#' @param loss `function(theta)` returning the scalar loss to MINIMISE (e.g. the
#'   distance between simulated and observed summary statistics).
#' @param ranges Named list of `c(lower, upper)` bounds, one per parameter.
#' @param n_init Number of initial random design points.
#' @param n_iter Number of Bayesian-optimisation iterations (loss evaluations after
#'   the initial design).
#' @param xi Expected-improvement exploration margin.
#' @param n_cand Candidate points scored by expected improvement each iteration.
#' @param seed Optional RNG seed.
#' @return a `gp_calibration`: the best `par` and `value`, the full design
#'   `X`/`y`, the best-so-far `trace`, and the final `gp` surrogate.
#' @references Jones DR, Schonlau M, Welch WJ (1998) J Glob Optim 13:455-492
#'   (efficient global optimisation / expected improvement).
#' @seealso [gpEmulator()], [abcSMC()], [personalizeCardio()]
#' @export
#' @examples
#' # recover the minimum of a 2-D loss in a few evaluations
#' loss <- function(p) (p[1] - 0.3)^2 + (p[2] - 0.7)^2
#' cal <- gpCalibrate(loss, list(x = c(0, 1), y = c(0, 1)), n_init = 8, n_iter = 15, seed = 1)
#' round(cal$par, 2)
gpCalibrate <- function(loss, ranges, n_init = 12, n_iter = 20, xi = 0.01,
                        n_cand = 500, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  ranges <- lapply(ranges, as.numeric); p <- length(ranges); nm <- names(ranges)
  if (is.null(nm)) nm <- paste0("p", seq_len(p))
  draw <- function(n) matrix(vapply(ranges, function(r) stats::runif(n, r[1], r[2]), numeric(n)),
                             n, p, dimnames = list(NULL, nm))
  X <- draw(n_init); y <- apply(X, 1, function(r) loss(stats::setNames(r, nm)))
  trace <- numeric(n_iter)
  for (it in seq_len(n_iter)) {
    gp <- gpEmulator(X, y, optimize = TRUE)
    Xc <- draw(n_cand); pr <- stats::predict(gp, Xc); fbest <- min(y)
    s <- pmax(pr$sd, 1e-9); z <- (fbest - pr$mean - xi) / s
    ei <- (fbest - pr$mean - xi) * stats::pnorm(z) + s * stats::dnorm(z)   # expected improvement
    ei[pr$sd < 1e-12] <- 0
    xnew <- Xc[which.max(ei), , drop = FALSE]
    ynew <- loss(stats::setNames(as.numeric(xnew), nm))
    X <- rbind(X, xnew); y <- c(y, ynew); trace[it] <- min(y)
  }
  bi <- which.min(y)
  structure(list(par = stats::setNames(X[bi, ], nm), value = y[bi], X = X, y = y,
                 trace = trace, gp = gpEmulator(X, y, optimize = TRUE), pnames = nm),
            class = "gp_calibration")
}

#' @export
print.gp_calibration <- function(x, ...) {
  cat(sprintf("GP calibration -- %d evaluations, best loss %.4g\n", length(x$y), x$value))
  cat(sprintf("  %s\n", paste(sprintf("%s = %.3g", x$pnames, x$par), collapse = ", ")))
  invisible(x)
}
