# Layer 3 -- data assimilation: fuse the mechanistic twin with real measurements.
#
# The unscented Kalman filter is a general nonlinear state-space estimator: given a
# process model f (how the state evolves) and a measurement model h (what the
# sensor sees), it tracks the state -- and, by AUGMENTING the state with unknown
# parameters, personalises the twin to a recorded sensor stream. This is what
# turns a generic simulator into a subject-specific digital twin.

# matrix square root via Cholesky with jitter fallback (keeps sigma points valid)
.ukf_sqrt <- function(M) {
  M <- (M + t(M)) / 2
  out <- tryCatch(t(chol(M)), error = function(e) NULL)
  if (!is.null(out)) return(out)
  e <- eigen(M, symmetric = TRUE); d <- pmax(e$values, 1e-10)
  e$vectors %*% diag(sqrt(d), length(d)) %*% t(e$vectors)
}

#' Unscented Kalman filter (general nonlinear state estimator)
#'
#' A sigma-point (Julier-Uhlmann) Kalman filter for a nonlinear state-space model
#' `x_k = f(x_{k-1}, k) + w`, `y_k = h(x_k, k) + v`. Used here to assimilate a
#' mechanistic twin with sensor data and to estimate augmented parameters, but it
#' is a stand-alone engine for any such model.
#'
#' @param f Process model `function(x, k)` returning the next state.
#' @param h Measurement model `function(x, k)` returning the predicted measurement.
#' @param y Measurement matrix (`n_steps x m`).
#' @param Q Process-noise covariance (state x state).
#' @param R Measurement-noise covariance (m x m).
#' @param x0,P0 Initial state mean and covariance.
#' @param alpha,beta,kappa Sigma-point spread parameters (defaults 1e-3, 2, 0).
#' @return a `ukf_result`: `state` (`n_steps x L` filtered means), `cov` (final
#'   covariance), `innovation` (`n_steps x m`).
#' @references Julier SJ, Uhlmann JK (1997) SPIE 3068:182-193; Wan & van der Merwe
#'   (2000).
#' @seealso [personalizeTwin()]
#' @export
#' @examples
#' # track a noisy sine with a random-walk state
#' set.seed(1); y <- matrix(sin(seq(0, 6, length.out = 60)) + rnorm(60, 0, .1), 60, 1)
#' fit <- unscentedKalmanFilter(function(x, k) x, function(x, k) x, y,
#'          Q = matrix(0.05), R = matrix(0.01), x0 = 0, P0 = matrix(1))
#' fit$state[60]
unscentedKalmanFilter <- function(f, h, y, Q, R, x0, P0,
                                  alpha = 1e-3, beta = 2, kappa = 0) {
  y <- as.matrix(y); nsteps <- nrow(y); m <- ncol(y)
  L <- length(x0); lambda <- alpha^2 * (L + kappa) - L
  Wm <- Wc <- numeric(2 * L + 1)
  Wm[1] <- lambda / (L + lambda); Wc[1] <- Wm[1] + (1 - alpha^2 + beta)
  Wm[-1] <- Wc[-1] <- 1 / (2 * (L + lambda))
  sig <- function(x, P) { S <- .ukf_sqrt((L + lambda) * P); cbind(x, x + S, x - S) }
  X <- matrix(0, nsteps, L); innov <- matrix(0, nsteps, m)
  x <- as.numeric(x0); P <- P0
  for (k in seq_len(nsteps)) {
    sp <- sig(x, P)
    spf <- matrix(vapply(seq_len(2 * L + 1), function(j) f(sp[, j], k), numeric(L)), nrow = L)
    xp <- as.numeric(spf %*% Wm)
    Pp <- Q; for (j in seq_len(2 * L + 1)) { d <- spf[, j] - xp; Pp <- Pp + Wc[j] * tcrossprod(d) }
    sp2 <- sig(xp, Pp)
    spy <- vapply(seq_len(2 * L + 1), function(j) h(sp2[, j], k), numeric(m))
    spy <- matrix(spy, nrow = m)
    yp <- as.numeric(spy %*% Wm)
    Pyy <- R; Pxy <- matrix(0, L, m)
    for (j in seq_len(2 * L + 1)) {
      dy <- spy[, j] - yp; dx <- sp2[, j] - xp
      Pyy <- Pyy + Wc[j] * tcrossprod(dy); Pxy <- Pxy + Wc[j] * outer(dx, dy)
    }
    K <- Pxy %*% solve(Pyy)
    resid <- y[k, ] - yp
    x <- xp + as.numeric(K %*% resid); P <- Pp - K %*% Pyy %*% t(K)
    X[k, ] <- x; innov[k, ] <- resid
  }
  structure(list(state = X, cov = P, innovation = innov, L = L), class = "ukf_result")
}

#' @export
print.ukf_result <- function(x, ...) {
  cat(sprintf("UKF -- %d steps, %d-D state, RMS innovation %.3g\n",
              nrow(x$state), x$L, sqrt(mean(x$innovation^2))))
  invisible(x)
}

#' Personalise a movement twin to recorded IMU data
#'
#' Estimates a subject's twin parameters by assimilating their IMU stream with an
#' unscented Kalman filter over an augmented state `[theta, omega, log(params)]`.
#' Parameters are estimated in log-space (keeping them positive) and the fitted
#' values are read from the converged filter.
#'
#' @param imu A data frame with `gyro, ax, ay` (as from [imuMeasure()]), sampled
#'   at `dt`.
#' @param prior A `movement_twin` giving the fixed parameters and the initial
#'   guess for the estimated ones.
#' @param estimate Character vector of parameters to personalise (any of
#'   `"strength"`, `"damping"`, `"stiffness"`, `"inertia"`).
#' @param method `"optim"` (default) fits the whole record by minimising the
#'   simulated-vs-observed IMU error (a batch/variational fit, robust to strong
#'   non-linearity); `"ukf"` runs the recursive unscented Kalman filter (online,
#'   best in the mildly non-linear regime).
#' @param dt Sample period (s).
#' @param param_rw Random-walk SD (log-space) for the estimated parameters
#'   (their assumed drift; default 0.02).
#' @param gyro_var,accel_var Measurement-noise variances for the gyro and accel
#'   channels.
#' @param burn_in Fraction of the record discarded before averaging the converged
#'   estimate (default 0.5).
#' @param n_starts Number of coarse-grid starting points scanned before the
#'   `"optim"` fit (multistart). The objective is evaluated on a multiplicative
#'   grid around the prior (one cheap simulation each) and the fit is started from
#'   the best point, making multi-parameter fits robust to spurious co-estimation
#'   minima at negligible cost; `1` reproduces a single Nelder-Mead start from the
#'   prior.
#' @return a `personalized_twin`: `twin` (the fitted twin), `estimates` (named
#'   fitted parameters), the co-estimated `gyro_bias`, and -- for `method =
#'   "ukf"` -- the per-step `trajectory` and the `ukf` result, or -- for `method
#'   = "optim"` -- the fitted `theta0` and final `objective`.
#' @seealso [unscentedKalmanFilter()], [extendedKalmanFilter()], [particleFilter()],
#'   [abcCalibration()], [insilicoIntervention()]
#' @export
#' @examples
#' \donttest{
#' truth <- limbTwin(strength = 7, damping = 0.5)
#' sim <- simulateTwin(truth, duration = 6, dt = 0.01)
#' imu <- imuMeasure(sim, imuSensor(gyro_noise = 0.03, accel_noise = 0.1), seed = 1)
#' fit <- personalizeTwin(imu, limbTwin(strength = 4, damping = 0.2),
#'                        estimate = c("strength", "damping"), dt = 0.01)
#' fit$estimates
#' }
personalizeTwin <- function(imu, prior, estimate = "strength", method = c("optim", "ukf"),
                            dt = 0.005, param_rw = 0.02, gyro_var = 0.03^2,
                            accel_var = 0.15^2, burn_in = 0.5, n_starts = 25L) {
  stopifnot(inherits(prior, "movement_twin"))
  method <- match.arg(method)
  estimate <- match.arg(estimate, c("strength", "damping", "stiffness", "inertia"),
                        several.ok = TRUE)
  np <- length(estimate); p0 <- unlist(prior[estimate])
  mk_est <- function(z) { tw <- prior; tw[estimate] <- as.list(exp(z)); tw }
  if (method == "optim") {
    Yg <- imu$gyro; Ya <- imu$ax; Yb <- imu$ay; n <- length(Yg)
    dur <- (n - 1) * dt; keep <- floor(0.2 * n):n     # skip the initial transient
    obj <- function(par) {
      tw <- mk_est(par[seq_len(np)]); theta0 <- par[np + 1]; gbias <- par[np + 2]
      sim <- simulateTwin(tw, dur, dt, theta0 = theta0, omega0 = Yg[1] - gbias)
      pred <- imuMeasure(sim, noise = FALSE)
      sum(((pred$gyro[keep] + gbias) - Yg[keep])^2) / gyro_var +
        sum((pred$ax[keep] - Ya[keep])^2 + (pred$ay[keep] - Yb[keep])^2) / accel_var
    }
    lp0 <- log(p0)
    if (n_starts <= 1L) {
      opt <- stats::optim(c(lp0, 0, 0), obj, method = "Nelder-Mead",
                          control = list(maxit = 800, reltol = 1e-8))
    } else {
      # Multistart by a cheap coarse grid over the estimated parameters (a
      # multiplicative grid around the prior). The fit objective is evaluated at
      # each grid point -- one simulation each, far cheaper than an optim -- and a
      # single Nelder-Mead is then run from the best-fitting grid point. This finds
      # the correct basin when a lone start would fall into a spurious co-estimation
      # minimum, at negligible extra cost. Deterministic, so still reproducible.
      levels <- log(c(0.3, 0.55, 1, 1.8, 3.3))
      grid <- as.matrix(expand.grid(rep(list(levels), np)))
      if (nrow(grid) > n_starts)
        grid <- grid[round(seq(1, nrow(grid), length.out = n_starts)), , drop = FALSE]
      go <- apply(grid, 1, function(g) obj(c(lp0 + g, 0, 0)))
      opt <- stats::optim(c(lp0 + grid[which.min(go), ], 0, 0), obj, method = "Nelder-Mead",
                          control = list(maxit = 500, reltol = 1e-8))  # good grid start -> fewer iters
    }
    est <- exp(opt$par[seq_len(np)]); names(est) <- estimate
    return(structure(list(twin = mk_est(log(est)), estimates = est,
                          gyro_bias = opt$par[np + 2], theta0 = opt$par[np + 1],
                          objective = opt$value, method = "optim", estimated = estimate),
                     class = "personalized_twin"))
  }
  # a gyroscope-bias nuisance state is ALWAYS co-estimated so a constant bias /
  # slow drift in the IMU does not corrupt the parameter estimates.
  Y <- as.matrix(imu[, c("gyro", "ax", "ay")])
  bi <- 2 + np + 1L                                    # index of the gyro-bias state
  mk_twin <- function(z) { tw <- prior; tw[estimate] <- as.list(exp(z[seq_len(np)])); tw }
  f <- function(x, k) {                                # RK4 one step of the dynamics
    tw <- mk_twin(x[-(1:2)]); s <- x[1:2]; t <- (k - 1) * dt
    k1 <- .twin_deriv(t, s, tw); k2 <- .twin_deriv(t + dt/2, s + dt/2*k1, tw)
    k3 <- .twin_deriv(t + dt/2, s + dt/2*k2, tw); k4 <- .twin_deriv(t + dt, s + dt*k3, tw)
    c(s + dt/6 * (k1 + 2*k2 + 2*k3 + k4), x[-(1:2)])   # params + bias: random walk
  }
  h <- function(x, k) {                                # IMU forward model + gyro bias
    tw <- mk_twin(x[-(1:2)]); t <- (k - 1) * dt
    alpha <- .twin_deriv(t, x[1:2], tw)[2]
    m <- .imu_forward(x[1], x[2], alpha, prior$sensor_r, prior$g)
    m["gyro"] <- m["gyro"] + x[bi]; m
  }
  x0 <- c(0, Y[1, 1], log(p0), 0)
  P0 <- diag(c(0.2, 0.5, rep(0.5, np), 0.1))
  Q  <- diag(c(1e-6, 1e-5, rep(param_rw^2, np), 1e-7)) # bias drifts slowly
  R  <- diag(c(gyro_var, accel_var, accel_var))
  uk <- unscentedKalmanFilter(f, h, Y, Q, R, x0, P0)
  traj <- exp(uk$state[, 2 + seq_len(np), drop = FALSE]); colnames(traj) <- estimate
  keep <- seq(max(1, floor(burn_in * nrow(traj))), nrow(traj))
  est <- colMeans(traj[keep, , drop = FALSE])
  structure(list(twin = mk_twin(log(est)), estimates = est, trajectory = traj,
                 gyro_bias = mean(uk$state[keep, bi]), ukf = uk, estimated = estimate),
            class = "personalized_twin")
}

#' @export
print.personalized_twin <- function(x, ...) {
  cat("Personalised movement twin -- estimated:",
      paste(sprintf("%s=%.3f", names(x$estimates), x$estimates), collapse = ", "), "\n")
  invisible(x)
}
