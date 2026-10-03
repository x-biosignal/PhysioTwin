# Layer 3/4 -- whole-record Bayesian fitting of the closed-loop cardiovascular twin.
#
# personalizeCardio() fits the closed loop to a few HRV/BPV SUMMARIES and the
# parameters are only weakly identified (a ridge collapses gains that trade off).
# The WHOLE RR record carries far more: the Mayer (low-frequency) rhythm and the
# respiratory (high-frequency) rhythm sit at different frequencies, so their gains
# are separately identifiable from the time series even though the scalar LF/HF
# powers confound them. cardioStateSpace() casts the closed loop as a tractable
# linear-Gaussian state-space model -- a damped Mayer oscillator (the delayed
# baroreflex resonance) plus a respiratory drive -- and personalizeCardioWaveform()
# fits its parameters jointly from the whole record by particle MCMC. Base R.

#' Closed-loop cardiovascular state-space model (for whole-record fitting)
#'
#' A linear-Gaussian state-space representation of the RR-interval series: a
#' damped second-order Mayer oscillator (the delayed-baroreflex resonance) drives
#' the heart period through a gain `g`, and respiration adds a high-frequency
#' (respiratory sinus arrhythmia) term of amplitude `A_rsa`. Suitable as the
#' `build_model` argument of [particleMCMC()] / [particleGibbs()].
#'
#' @param g Low-frequency gain onto the RR interval (ms per unit state).
#' @param rho Mayer-oscillator damping (0..1; nearer 1 = sharper low-frequency peak;
#'   used when `lf_model = "oscillator"`).
#' @param omega_m Mayer angular frequency (radians per beat; `"oscillator"` only).
#' @param A_rsa Respiratory (RSA) amplitude on the RR interval (ms).
#' @param resp Respiratory regressor sampled at the beats (length = number of beats).
#' @param rr0 Mean RR interval (ms).
#' @param q_mayer,r_obs Low-frequency process-noise and RR observation-noise variances.
#' @param lf_model Low-frequency model: `"oscillator"` (a resonant, second-order
#'   damped Mayer oscillator, for records with a distinct Mayer peak) or `"ar1"` (a
#'   first-order autoregressive process giving BROADBAND low-frequency power without a
#'   peak, for records whose low-frequency variability is not a sharp rhythm).
#' @param phi First-order autoregressive coefficient (0..1; used when
#'   `lf_model = "ar1"`); nearer 1 = more low-frequency power.
#' @return a state-space model list `f`, `h` (`function(x, k)`), `Q`, `R`, `x0`, `P0`.
#' @seealso [personalizeCardioWaveform()], [particleMCMC()], [cardioRespiratory()]
#' @export
#' @examples
#' resp <- sin(2 * pi * 0.25 * cumsum(rep(0.85, 200)))
#' m <- cardioStateSpace(g = 8, rho = 0.9, omega_m = 0.53, A_rsa = 20, resp = resp)
#' m$h(c(1, 0), 1)
cardioStateSpace <- function(g = 8, rho = 0.9, omega_m = 0.53, A_rsa = 20,
                             resp, rr0 = 850, q_mayer = 1, r_obs = 9,
                             lf_model = c("oscillator", "ar1"), phi = 0.9) {
  lf_model <- match.arg(lf_model)
  if (lf_model == "ar1") {                                     # broadband LF: first-order AR, no peak
    s2 <- q_mayer / (1 - phi^2)
    return(list(f = function(x, k) phi * x[1],
                h = function(x, k) rr0 + g * x[1] + A_rsa * resp[k],
                Q = matrix(q_mayer), R = matrix(r_obs), x0 = 0, P0 = matrix(s2),
                lf_model = "ar1"))
  }
  A <- matrix(c(2 * rho * cos(omega_m), 1, -rho^2, 0), 2, 2)   # damped-oscillator transition
  s2 <- q_mayer / (1 - rho^2)                                  # stationary state variance
  list(f = function(x, k) as.numeric(A %*% x),
       h = function(x, k) rr0 + g * x[1] + A_rsa * resp[k],
       Q = matrix(c(q_mayer, 0, 0, 1e-8), 2, 2), R = matrix(r_obs),
       x0 = c(0, 0), P0 = diag(c(s2, s2), 2), lf_model = "oscillator")
}

# Exact Kalman marginal log-likelihood of an affine-Gaussian model
# y_k = offset_k + H x_k + v,  x_k = A x_{k-1} + w. Linear-Gaussian, so exact and
# fast -- no particles needed (particle MCMC is for the intractable, nonlinear case).
.kf_loglik <- function(A, H, offset, Q, R, x0, P0, y) {
  Tn <- nrow(y); m <- ncol(y); x <- x0; P <- P0; ll <- 0
  if (m == 1) {                                         # scalar-observation fast path (RR series)
    r <- R[1, 1]
    for (k in seq_len(Tn)) {
      xp <- A %*% x; Pp <- A %*% P %*% t(A) + Q; Hp <- H %*% Pp
      S <- as.numeric(Hp %*% t(H)) + r; v <- y[k] - offset[k] - as.numeric(H %*% xp)
      ll <- ll - 0.5 * (v * v / S + log(2 * pi * S))
      K <- t(Hp) / S; x <- xp + K * v; P <- Pp - K %*% Hp
    }
    return(ll)
  }
  for (k in seq_len(Tn)) {
    xp <- A %*% x; Pp <- A %*% P %*% t(A) + Q
    v <- y[k, ] - offset[k] - H %*% xp
    S <- H %*% Pp %*% t(H) + R; Sinv <- solve(S)
    ll <- ll - 0.5 * (as.numeric(crossprod(v, Sinv %*% v)) +
                        as.numeric(determinant(S, logarithm = TRUE)$modulus) + m * log(2 * pi))
    K <- Pp %*% t(H) %*% Sinv; x <- xp + K %*% v; P <- Pp - K %*% H %*% Pp
  }
  ll
}

#' Personalise the closed-loop cardiovascular twin from the whole RR record
#'
#' Fits the [cardioStateSpace()] parameters jointly to a whole RR-interval record.
#' Because the closed-loop model is linear-Gaussian its likelihood is exact (a
#' Kalman filter), so the posterior is sampled by [metropolis()] -- the
#' whole-waveform counterpart to the summary-based [personalizeCardio()]. The Mayer
#' gain/damping and the respiratory amplitude are identifiable here because the
#' record separates the low- and high-frequency rhythms, resolving the ridge that
#' confounds the summary fit. (For an intractable, nonlinear model use
#' [particleMCMC()] / [particleGibbs()] with [cardioStateSpace()]'s `f`/`h`.)
#'
#' @param rr RR-interval series (ms).
#' @param resp_freq Respiratory frequency (Hz) for the RSA regressor.
#' @param estimate Parameters to fit (any of `"g"`, `"rho"`, `"A_rsa"`, `"phi"`);
#'   `NULL` picks the default for the `lf_model` (`g`/`rho`/`A_rsa` for the
#'   oscillator, `g`/`phi`/`A_rsa` for `"ar1"`).
#' @param lf_model Low-frequency model (see [cardioStateSpace()]): `"oscillator"`
#'   (a resonant Mayer peak) or `"ar1"` (broadband low-frequency power, for records
#'   without a distinct Mayer rhythm, e.g. a regular healthy sinus rhythm).
#' @param fixed Named list of held-fixed [cardioStateSpace()] parameters
#'   (`omega_m`, `rr0`, `q_mayer`, `r_obs`, and any low-frequency parameter not fitted).
#' @param n_iter Metropolis iterations.
#' @param seed Optional RNG seed.
#' @return a `cardio_waveform_fit`: the `posterior` (an `mcmc_result`), the fitted
#'   `estimates`, and the `estimate`d names.
#' @seealso [cardioStateSpace()], [metropolis()], [particleMCMC()], [personalizeCardio()]
#' @export
#' @examples
#' \donttest{
#' truth <- cardioStateSpace(g = 8, rho = 0.9, A_rsa = 20,
#'                           resp = sin(2 * pi * 0.25 * cumsum(rep(0.85, 250))))
#' set.seed(1); x <- c(0, 0); rr <- numeric(250)
#' for (k in seq_len(250)) {
#'   x <- truth$f(x, k) + c(rnorm(1), 0); rr[k] <- truth$h(x, k) + rnorm(1, 0, 3)
#' }
#' fit <- personalizeCardioWaveform(rr, resp_freq = 0.25, estimate = c("g", "A_rsa"),
#'                                  n_iter = 3000, seed = 1)
#' fit$estimates
#' }
personalizeCardioWaveform <- function(rr, resp_freq = 0.25, estimate = NULL,
                                      lf_model = c("oscillator", "ar1"),
                                      fixed = list(omega_m = 0.53, rr0 = NULL,
                                                   q_mayer = 1, r_obs = 9),
                                      n_iter = 4000, seed = NULL) {
  lf_model <- match.arg(lf_model); rr <- as.numeric(rr); Tn <- length(rr)
  if (is.null(estimate)) estimate <- if (lf_model == "ar1") c("g", "phi", "A_rsa") else c("g", "rho", "A_rsa")
  estimate <- match.arg(estimate, c("g", "rho", "A_rsa", "phi"), several.ok = TRUE)
  if (is.null(fixed$rr0)) fixed$rr0 <- mean(rr)
  Y <- matrix(rr, Tn, 1)
  # RSA regressor: respiration phase at the beats, on the mean-RR grid (avoids the
  # phase drift that a variable-interval grid would introduce over a long record)
  resp <- sin(2 * pi * resp_freq * (seq_len(Tn) - 1) * mean(rr) / 1000)
  bounds <- list(g = c(0.5, 300), rho = c(0.3, 0.995), A_rsa = c(0.1, 200), phi = c(0.3, 0.995))
  logpost <- function(theta) {                                  # exact Kalman likelihood + uniform prior
    if (!all(vapply(seq_along(estimate), function(j)
        theta[j] > bounds[[estimate[j]]][1] && theta[j] < bounds[[estimate[j]]][2], logical(1))))
      return(-Inf)
    p <- list(g = 8, rho = 0.9, A_rsa = 20, phi = 0.9)
    for (j in seq_along(estimate)) p[[estimate[j]]] <- theta[j]
    if (lf_model == "ar1") {                                    # broadband first-order-AR low frequency
      s2 <- fixed$q_mayer / (1 - p$phi^2)
      .kf_loglik(matrix(p$phi, 1, 1), matrix(p$g, 1, 1), fixed$rr0 + p$A_rsa * resp,
                 matrix(fixed$q_mayer), matrix(fixed$r_obs), 0, matrix(s2), Y)
    } else {                                                    # resonant Mayer oscillator
      A <- matrix(c(2 * p$rho * cos(fixed$omega_m), 1, -p$rho^2, 0), 2, 2)
      s2 <- fixed$q_mayer / (1 - p$rho^2)
      .kf_loglik(A, matrix(c(p$g, 0), 1, 2), fixed$rr0 + p$A_rsa * resp,
                 matrix(c(fixed$q_mayer, 0, 0, 1e-8), 2, 2), matrix(fixed$r_obs),
                 c(0, 0), diag(c(s2, s2), 2), Y)
    }
  }
  init <- stats::setNames(vapply(estimate, function(p)
    switch(p, g = 5, rho = 0.85, phi = 0.85, A_rsa = 10), numeric(1)), estimate)
  psd <- vapply(estimate, function(p) 0.05 * diff(bounds[[p]]), numeric(1))
  post <- metropolis(logpost, init, n_iter = n_iter, proposal_sd = psd, seed = seed)
  structure(list(posterior = post, estimates = post$mean, estimate = estimate, lf_model = lf_model),
            class = "cardio_waveform_fit")
}

#' @export
print.cardio_waveform_fit <- function(x, ...) {
  cat("Closed-loop whole-record fit (particle MCMC):\n")
  for (j in seq_along(x$estimate))
    cat(sprintf("  %-6s %.3f  [%.3f, %.3f]\n", x$estimate[j], x$posterior$mean[j],
                x$posterior$quantiles[1, j], x$posterior$quantiles[3, j]))
  invisible(x)
}
