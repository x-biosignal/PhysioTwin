# Layer 3 -- data assimilation: alternative nonlinear state-space estimators.
#
# The unscented Kalman filter (assimilate-ukf.R) is one way to fuse a mechanistic
# twin with measurements; different problems want different engines. All three
# here share the UKF's interface -- a process model f(x, k), a measurement model
# h(x, k) and measurements y -- and return the same shape, so they are drop-in
# alternatives:
#   * extendedKalmanFilter()  -- first-order (Jacobian) linearisation; cheapest,
#                                exact for a linear-Gaussian model. Jacobians are
#                                computed by finite differences unless supplied.
#   * particleFilter()        -- a bootstrap sequential-importance-resampling
#                                filter; a sample-based posterior that handles
#                                strong nonlinearity and non-Gaussian noise.
#   * ensembleKalmanFilter()  -- a stochastic ensemble Kalman filter; the
#                                ensemble carries the covariance, scaling to
#                                high-dimensional states.
# Dependency-free base R. Reuses .ukf_sqrt() from assimilate-ukf.R.

# Central-difference Jacobian of a vector field fn(x, k) at x (rows = outputs).
.num_jac <- function(fn, x, k) {
  f0 <- fn(x, k); m <- length(f0); L <- length(x); J <- matrix(0, m, L)
  for (i in seq_len(L)) {
    hi <- 1e-6 * max(1, abs(x[i]))
    xp <- x; xp[i] <- x[i] + hi; xm <- x; xm[i] <- x[i] - hi
    J[, i] <- (fn(xp, k) - fn(xm, k)) / (2 * hi)
  }
  J
}

# Draw n samples from N(mu, Sigma) (rows = samples).
.rmvn <- function(n, mu, Sigma) {
  L <- length(mu); S <- .ukf_sqrt(Sigma)
  Z <- matrix(stats::rnorm(n * L), n, L)
  sweep(Z %*% t(S), 2, mu, "+")
}

#' Extended Kalman filter (Jacobian-linearised state estimator)
#'
#' A first-order extended Kalman filter for `x_k = f(x_{k-1}, k) + w`,
#' `y_k = h(x_k, k) + v`. The process and measurement Jacobians are obtained by
#' central finite differences unless analytic versions are supplied. Same
#' interface and return shape as [unscentedKalmanFilter()]; for a linear-Gaussian
#' model it reproduces the exact Kalman filter.
#'
#' @param f,h Process and measurement models `function(x, k)`.
#' @param y Measurement matrix (`n_steps x m`).
#' @param Q,R Process- and measurement-noise covariances.
#' @param x0,P0 Initial state mean and covariance.
#' @param Fjac,Hjac Optional analytic Jacobians `function(x, k)` (rows = outputs);
#'   finite differences are used when `NULL`.
#' @return an `ekf_result`: `state` (`n_steps x L`), `cov` (final), `innovation`.
#' @references Jazwinski AH (1970) Stochastic Processes and Filtering Theory.
#' @seealso [unscentedKalmanFilter()], [particleFilter()], [ensembleKalmanFilter()]
#' @export
#' @examples
#' set.seed(1); y <- matrix(cumsum(rnorm(40, 0, .1)) + rnorm(40, 0, .1), 40, 1)
#' fit <- extendedKalmanFilter(function(x, k) x, function(x, k) x, y,
#'          Q = matrix(0.05), R = matrix(0.01), x0 = 0, P0 = matrix(1))
#' fit$state[40]
extendedKalmanFilter <- function(f, h, y, Q, R, x0, P0, Fjac = NULL, Hjac = NULL) {
  y <- as.matrix(y); nsteps <- nrow(y); m <- ncol(y); L <- length(x0)
  if (is.null(Fjac)) Fjac <- function(x, k) .num_jac(f, x, k)
  if (is.null(Hjac)) Hjac <- function(x, k) .num_jac(h, x, k)
  X <- matrix(0, nsteps, L); innov <- matrix(0, nsteps, m)
  x <- as.numeric(x0); P <- P0; Id <- diag(L)
  for (k in seq_len(nsteps)) {
    Fk <- Fjac(x, k); xp <- as.numeric(f(x, k)); Pp <- Fk %*% P %*% t(Fk) + Q
    Hk <- Hjac(xp, k); yp <- as.numeric(h(xp, k))
    S <- Hk %*% Pp %*% t(Hk) + R
    K <- Pp %*% t(Hk) %*% solve(S)
    resid <- y[k, ] - yp
    x <- xp + as.numeric(K %*% resid); P <- (Id - K %*% Hk) %*% Pp
    X[k, ] <- x; innov[k, ] <- resid
  }
  structure(list(state = X, cov = P, innovation = innov, L = L), class = "ekf_result")
}

#' @export
print.ekf_result <- function(x, ...) {
  cat(sprintf("EKF -- %d steps, %d-D state, RMS innovation %.3g\n",
              nrow(x$state), x$L, sqrt(mean(x$innovation^2)))); invisible(x)
}

# systematic resampling of a normalised weight vector -> parent indices
.sys_resample <- function(w) {
  N <- length(w); pos <- (stats::runif(1) + 0:(N - 1)) / N
  cw <- cumsum(w); idx <- integer(N); i <- 1L; j <- 1L
  while (i <= N) { if (pos[i] <= cw[j]) { idx[i] <- j; i <- i + 1L } else j <- j + 1L }
  idx
}

#' Bootstrap particle filter (sequential Monte-Carlo state estimator)
#'
#' A sequential-importance-resampling (bootstrap) particle filter for
#' `x_k = f(x_{k-1}, k) + w`, `y_k = h(x_k, k) + v`. Particles are propagated
#' through `f` with process noise `Q`, weighted by the Gaussian measurement
#' likelihood, and resampled (systematic) when the effective sample size falls
#' below `resample_thresh` of the population. Handles strong nonlinearity and
#' multi-modal posteriors that a Kalman filter cannot. Same interface as
#' [unscentedKalmanFilter()].
#'
#' @param f,h Process and measurement models `function(x, k)`.
#' @param y Measurement matrix (`n_steps x m`).
#' @param Q,R Process- and measurement-noise covariances.
#' @param x0,P0 Initial state mean and covariance.
#' @param n_particles Number of particles.
#' @param resample_thresh Resample when effective sample size `< thresh * N`.
#' @param seed Optional RNG seed.
#' @return a `pf_result`: `state` (`n_steps x L` posterior means), `cov` (final
#'   weighted covariance), `innovation`, and `ess` (effective sample size per step).
#' @references Gordon NJ, Salmond DJ, Smith AFM (1993) IEE Proc F 140:107-113.
#' @seealso [unscentedKalmanFilter()], [extendedKalmanFilter()], [ensembleKalmanFilter()]
#' @export
#' @examples
#' set.seed(1); y <- matrix(sin(seq(0, 6, length.out = 50)) + rnorm(50, 0, .1), 50, 1)
#' fit <- particleFilter(function(x, k) x, function(x, k) x, y,
#'          Q = matrix(0.05), R = matrix(0.01), x0 = 0, P0 = matrix(1), n_particles = 500)
#' fit$state[50]
particleFilter <- function(f, h, y, Q, R, x0, P0, n_particles = 1000,
                           resample_thresh = 0.5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  y <- as.matrix(y); nsteps <- nrow(y); m <- ncol(y); L <- length(x0)
  Rinv <- solve(R)
  parts <- t(.rmvn(n_particles, as.numeric(x0), P0))         # L x N
  w <- rep(1 / n_particles, n_particles)
  X <- matrix(0, nsteps, L); innov <- matrix(0, nsteps, m); ess <- numeric(nsteps)
  for (k in seq_len(nsteps)) {
    prop <- vapply(seq_len(n_particles), function(i) f(parts[, i], k), numeric(L))
    prop <- matrix(prop, nrow = L) + t(.rmvn(n_particles, rep(0, L), Q))
    hy <- matrix(vapply(seq_len(n_particles), function(i) h(prop[, i], k), numeric(m)), nrow = m)
    innov[k, ] <- y[k, ] - rowMeans(hy)                       # prior-mean prediction error
    d <- hy - y[k, ]
    ll <- -0.5 * colSums((Rinv %*% d) * d)
    w <- w * exp(ll - max(ll)); w <- w / sum(w)
    X[k, ] <- prop %*% w                                      # posterior mean
    ess[k] <- 1 / sum(w^2)
    parts <- prop
    if (ess[k] < resample_thresh * n_particles) {
      parts <- parts[, .sys_resample(w), drop = FALSE]; w <- rep(1 / n_particles, n_particles)
    }
  }
  dev <- parts - X[nsteps, ]
  P <- (dev %*% (t(dev) * w)) * (n_particles / (n_particles - 1))
  structure(list(state = X, cov = P, innovation = innov, ess = ess, L = L),
            class = "pf_result")
}

#' @export
print.pf_result <- function(x, ...) {
  cat(sprintf("Particle filter -- %d steps, %d-D state, mean ESS %.0f, RMS innovation %.3g\n",
              nrow(x$state), x$L, mean(x$ess), sqrt(mean(x$innovation^2)))); invisible(x)
}

#' Ensemble Kalman filter (Monte-Carlo covariance state estimator)
#'
#' A stochastic ensemble Kalman filter for `x_k = f(x_{k-1}, k) + w`,
#' `y_k = h(x_k, k) + v`. An ensemble of members carries the state distribution;
#' its sample cross- and innovation-covariances form the Kalman gain, and the
#' members are updated with perturbed observations. Scales to high-dimensional
#' states where a full covariance is impractical. Same interface as
#' [unscentedKalmanFilter()].
#'
#' @param f,h Process and measurement models `function(x, k)`.
#' @param y Measurement matrix (`n_steps x m`).
#' @param Q,R Process- and measurement-noise covariances.
#' @param x0,P0 Initial state mean and covariance.
#' @param n_ensemble Ensemble size.
#' @param seed Optional RNG seed.
#' @return an `enkf_result`: `state` (`n_steps x L` ensemble means), `cov` (final
#'   ensemble covariance), `innovation`.
#' @references Evensen G (1994) J Geophys Res 99:10143-10162; Burgers et al. (1998).
#' @seealso [unscentedKalmanFilter()], [extendedKalmanFilter()], [particleFilter()]
#' @export
#' @examples
#' set.seed(1); y <- matrix(sin(seq(0, 6, length.out = 50)) + rnorm(50, 0, .1), 50, 1)
#' fit <- ensembleKalmanFilter(function(x, k) x, function(x, k) x, y,
#'          Q = matrix(0.05), R = matrix(0.01), x0 = 0, P0 = matrix(1), n_ensemble = 80)
#' fit$state[50]
ensembleKalmanFilter <- function(f, h, y, Q, R, x0, P0, n_ensemble = 100, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  y <- as.matrix(y); nsteps <- nrow(y); m <- ncol(y); L <- length(x0); N <- n_ensemble
  E <- t(.rmvn(N, as.numeric(x0), P0))                        # L x N
  X <- matrix(0, nsteps, L); innov <- matrix(0, nsteps, m)
  for (k in seq_len(nsteps)) {
    E <- matrix(vapply(seq_len(N), function(i) f(E[, i], k), numeric(L)), nrow = L)
    E <- E + t(.rmvn(N, rep(0, L), Q))                        # forecast + process noise
    Hy <- matrix(vapply(seq_len(N), function(i) h(E[, i], k), numeric(m)), nrow = m)
    xbar <- rowMeans(E); ybar <- rowMeans(Hy)
    A <- E - xbar; B <- Hy - ybar
    Pxy <- (A %*% t(B)) / (N - 1); Pyy <- (B %*% t(B)) / (N - 1) + R
    K <- Pxy %*% solve(Pyy)
    obs <- y[k, ] + t(.rmvn(N, rep(0, m), R))                 # perturbed observations
    E <- E + K %*% (obs - Hy)
    X[k, ] <- rowMeans(E); innov[k, ] <- y[k, ] - ybar
  }
  dev <- E - rowMeans(E); P <- (dev %*% t(dev)) / (N - 1)
  structure(list(state = X, cov = P, innovation = innov, L = L), class = "enkf_result")
}

#' @export
print.enkf_result <- function(x, ...) {
  cat(sprintf("Ensemble KF -- %d steps, %d-D state, RMS innovation %.3g\n",
              nrow(x$state), x$L, sqrt(mean(x$innovation^2)))); invisible(x)
}
