# Layer 3 -- data assimilation: a Gaussian-process surrogate (emulator).
#
# A mechanistic twin can be expensive to run, which makes calibration, sensitivity
# analysis and optimisation slow. A Gaussian-process emulator learns a cheap
# probabilistic surrogate of the parameter -> output map from a modest design of
# simulator runs: it interpolates the training runs exactly and, crucially,
# reports its own predictive UNCERTAINTY -- large away from the training points,
# small near them -- so downstream methods know where the surrogate can be trusted
# and where the real simulator must be run. Squared-exponential kernel, marginal-
# likelihood hyperparameter fit, Cholesky solves. Dependency-free base R.

# squared Euclidean distances between rows of X and rows of Z
.sqdist <- function(X, Z) {
  outer(rowSums(X^2), rowSums(Z^2), "+") - 2 * X %*% t(Z)
}

# Cholesky with escalating jitter -- near-duplicate design points (e.g. a Bayesian
# optimiser sampling repeatedly near the optimum) make the kernel matrix singular
.safe_chol <- function(K) {
  n <- nrow(K); base <- 1e-10 * mean(diag(K)); jit <- base
  for (t in 0:8) {
    L <- tryCatch(chol(K + diag(jit, n)), error = function(e) NULL)
    if (!is.null(L)) return(L)
    jit <- max(jit * 10, 1e-12)
  }
  chol(K + diag(1e-6 * mean(diag(K)) + 1e-10, n))
}

# negative log marginal likelihood at log-hyperparameters (ell, sf, sn)
.gp_nlml <- function(par, D2, y) {
  ell <- exp(par[1]); sf2 <- exp(2 * par[2]); sn2 <- exp(2 * par[3]); n <- length(y)
  K <- sf2 * exp(-D2 / (2 * ell^2)) + diag(sn2 + 1e-8, n)
  L <- tryCatch(chol(K), error = function(e) NULL); if (is.null(L)) return(1e10)
  a <- backsolve(L, forwardsolve(t(L), y))
  0.5 * sum(y * a) + sum(log(diag(L))) + 0.5 * n * log(2 * pi)
}

#' Gaussian-process emulator (surrogate model)
#'
#' Fits a Gaussian-process regression surrogate with a squared-exponential kernel
#' to training inputs `X` and outputs `y` -- typically a design of simulator runs,
#' so the emulator stands in for an expensive twin. Hyperparameters (length-scale,
#' signal and noise SD) are fitted by maximising the marginal likelihood.
#' [predict.gpEmulator()] returns the posterior mean and standard deviation at new
#' inputs; the SD is small near training points and grows away from them.
#'
#' @param X Training inputs: a matrix (`n x d`) or length-`n` vector.
#' @param y Training outputs (length `n`).
#' @param optimize If `TRUE`, fit hyperparameters by marginal likelihood; otherwise
#'   use the initial values.
#' @param lengthscale,sigma_f,sigma_n Initial (or fixed) kernel length-scale,
#'   signal SD and noise SD; `NULL` picks data-driven defaults (median pairwise
#'   distance, `sd(y)`, a small nugget).
#' @return a `gpEmulator` object (with a [predict.gpEmulator()] method): the fitted
#'   `lengthscale`, `sigma_f`, `sigma_n` and the training data.
#' @references Rasmussen CE, Williams CKI (2006) Gaussian Processes for Machine Learning.
#' @seealso [predict.gpEmulator()], [sobolIndices()], [abcCalibration()]
#' @export
#' @examples
#' set.seed(1); x <- runif(20, 0, 6); y <- sin(x) + rnorm(20, 0, 0.05)
#' gp <- gpEmulator(x, y)
#' p <- predict(gp, seq(0, 6, length.out = 5)); round(p$mean, 2)
gpEmulator <- function(X, y, optimize = TRUE, lengthscale = NULL,
                       sigma_f = NULL, sigma_n = NULL) {
  X <- if (is.matrix(X)) X else matrix(X, ncol = 1)
  y <- as.numeric(y); n <- length(y); ybar <- mean(y); yc <- y - ybar
  D2 <- .sqdist(X, X)
  if (is.null(lengthscale)) { md <- sqrt(D2[upper.tri(D2)]); lengthscale <- stats::median(md[md > 0]) }
  if (is.null(sigma_f)) sigma_f <- stats::sd(yc); if (sigma_f < 1e-6) sigma_f <- 1
  if (is.null(sigma_n)) sigma_n <- max(1e-4, 0.05 * stats::sd(yc))
  par <- log(c(lengthscale, sigma_f, sigma_n))
  if (optimize) {
    opt <- tryCatch(stats::optim(par, .gp_nlml, D2 = D2, y = yc,
                                 method = "Nelder-Mead", control = list(maxit = 400)),
                    error = function(e) NULL)
    if (!is.null(opt)) par <- opt$par
  }
  ell <- exp(par[1]); sf <- exp(par[2]); sn <- exp(par[3])
  K <- sf^2 * exp(-D2 / (2 * ell^2)) + diag(sn^2 + 1e-8, n)
  L <- .safe_chol(K); alpha <- backsolve(L, forwardsolve(t(L), yc))
  structure(list(X = X, y = y, ybar = ybar, L = L, alpha = alpha,
                 lengthscale = ell, sigma_f = sf, sigma_n = sn),
            class = "gpEmulator")
}

#' Predict from a Gaussian-process emulator
#'
#' @param object a `gpEmulator` from [gpEmulator()].
#' @param newdata New inputs: a matrix (`m x d`) or a vector (for 1-D emulators).
#' @param ... Unused.
#' @return a list with `mean` (posterior mean) and `sd` (posterior standard
#'   deviation) at each new input.
#' @seealso [gpEmulator()]
#' @export
predict.gpEmulator <- function(object, newdata, ...) {
  Xn <- if (is.matrix(newdata)) newdata else matrix(newdata, ncol = ncol(object$X))
  ell <- object$lengthscale; sf <- object$sigma_f
  Ks <- sf^2 * exp(-.sqdist(Xn, object$X) / (2 * ell^2))            # m x n
  mu <- object$ybar + as.numeric(Ks %*% object$alpha)
  v <- forwardsolve(t(object$L), t(Ks))                            # n x m
  var <- sf^2 - colSums(v^2)
  list(mean = mu, sd = sqrt(pmax(0, var)))
}

#' @export
print.gpEmulator <- function(x, ...) {
  cat(sprintf("GP emulator -- %d training points, %d-D input, lengthscale %.3g, signal %.3g, noise %.3g\n",
              nrow(x$X), ncol(x$X), x$lengthscale, x$sigma_f, x$sigma_n)); invisible(x)
}
