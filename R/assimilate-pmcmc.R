# Layer 3 -- data assimilation: particle MCMC for whole-waveform Bayesian fitting.
#
# personalizeCardio() fits the closed loop to a few variability SUMMARIES, and the
# parameters are only weakly identified (a ridge). Fitting the WHOLE waveform
# instead uses all the information in the time series. For a state-space model
# whose likelihood is intractable, particle MCMC (particle marginal
# Metropolis-Hastings; Andrieu, Doucet & Holenstein 2010) is exact: a particle
# filter supplies an UNBIASED estimate of the marginal likelihood of the whole
# observation series given the parameters, and Metropolis-Hastings uses that
# estimate in the acceptance ratio -- targeting the true parameter posterior
# despite the noisy likelihood. Reuses the particle machinery of assimilate-filters.R.

# Bootstrap particle filter estimate of the log marginal likelihood log p(y_{1:T} | theta).
.pf_loglik <- function(f, h, y, Q, R, x0, P0, N) {
  Y <- as.matrix(y); Tn <- nrow(Y); m <- ncol(Y); L <- length(x0)
  Rinv <- solve(R); logdetR <- as.numeric(determinant(R, logarithm = TRUE)$modulus)
  X <- t(.rmvn(N, as.numeric(x0), P0))                      # L x N initial particles
  ll <- 0
  for (k in seq_len(Tn)) {
    X <- matrix(vapply(seq_len(N), function(i) f(X[, i], k), numeric(L)), L, N) +
           t(.rmvn(N, rep(0, L), Q))                        # propagate + process noise
    HX <- matrix(vapply(seq_len(N), function(i) h(X[, i], k), numeric(m)), m, N)
    v <- Y[k, ] - HX
    logw <- -0.5 * (colSums(v * (Rinv %*% v)) + m * log(2 * pi) + logdetR)
    mx <- max(logw); w <- exp(logw - mx)
    ll <- ll + mx + log(mean(w) + 1e-300)                   # unbiased incremental log-likelihood (floored)
    if (!is.finite(ll)) return(-.Machine$double.xmax)       # a finite floor lets the sampler escape
    sw <- sum(w); if (sw <= 0) return(-.Machine$double.xmax)
    X <- X[, .sys_resample(w / sw), drop = FALSE]           # systematic resampling
  }
  ll
}

#' Particle marginal Metropolis-Hastings (particle MCMC)
#'
#' Samples the parameter posterior of a state-space model from a WHOLE observation
#' series, when the likelihood is intractable. Each proposed parameter is scored by
#' a bootstrap particle filter's unbiased estimate of the marginal likelihood
#' `p(y_{1:T} | theta)`, which enters the Metropolis-Hastings acceptance ratio;
#' the resulting chain targets the exact parameter posterior (Andrieu-Doucet-
#' Holenstein). The proposal covariance adapts during burn-in.
#'
#' @param build_model `function(theta)` returning the state-space model for
#'   parameters `theta`: a list `f`, `h` (`function(x, k)`), `Q`, `R`, `x0`, `P0`.
#' @param y Observation matrix (`n_steps x m`) or vector.
#' @param prior_logd `function(theta)` returning the log prior density (`-Inf` off
#'   the support).
#' @param init Initial parameter vector (its names label the output).
#' @param n_particles Particles in the likelihood-estimating filter.
#' @param n_iter MCMC iterations.
#' @param proposal_sd Initial per-parameter proposal SD.
#' @param burn_in Fraction discarded (and over which the proposal adapts).
#' @param seed Optional RNG seed.
#' @return a `pmcmc_result`: posterior `samples`, `acceptance`, `mean`, `sd` and
#'   2.5/50/97.5% `quantiles`.
#' @references Andrieu C, Doucet A, Holenstein R (2010) J R Stat Soc B 72:269-342.
#' @seealso [particleFilter()], [metropolis()], [personalizeCardio()]
#' @export
#' @examples
#' \donttest{
#' # recover the AR coefficient of a linear-Gaussian state-space model from the series
#' set.seed(1); a <- 0.7; x <- numeric(60); for (k in 2:60) x[k] <- a * x[k-1] + rnorm(1, 0, 0.3)
#' y <- matrix(x + rnorm(60, 0, 0.3), 60, 1)
#' bm <- function(th) list(f = function(x, k) th[1] * x, h = function(x, k) x,
#'                         Q = matrix(0.09), R = matrix(0.09), x0 = 0, P0 = matrix(1))
#' fit <- particleMCMC(bm, y, function(th) dunif(th[1], -1, 1, log = TRUE),
#'                     init = c(a = 0.4), n_particles = 100, n_iter = 800, seed = 1)
#' fit$mean
#' }
particleMCMC <- function(build_model, y, prior_logd, init, n_particles = 200,
                         n_iter = 3000, proposal_sd = 0.1, burn_in = 0.3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  d <- length(init); nm <- names(init); if (is.null(nm)) nm <- paste0("p", seq_len(d))
  ll_at <- function(th) {
    m <- build_model(th); .pf_loglik(m$f, m$h, y, m$Q, m$R, m$x0, m$P0, n_particles)
  }
  theta <- as.numeric(init); ll <- ll_at(theta); lpri <- prior_logd(theta)
  chain <- matrix(0, n_iter, d); acc <- logical(n_iter)
  Lp <- diag(proposal_sd, d); bwin <- floor(burn_in * n_iter)
  for (i in seq_len(n_iter)) {
    thp <- theta + as.numeric(Lp %*% stats::rnorm(d)); lprip <- prior_logd(thp)
    if (is.finite(lprip)) {
      llp <- ll_at(thp)
      if (is.finite(llp) && log(stats::runif(1)) < (llp + lprip - ll - lpri)) {
        theta <- thp; ll <- llp; lpri <- lprip; acc[i] <- TRUE
      }
    }
    chain[i, ] <- theta
    if (i >= 100 && i <= bwin && i %% 50 == 0) {            # adapt proposal (burn-in only)
      C <- stats::cov(chain[max(1, i - 400):i, , drop = FALSE]) * (2.4^2 / d) + diag(1e-8, d)
      Lp <- tryCatch(t(chol(C)), error = function(e) Lp)
    }
  }
  keep <- (bwin + 1):n_iter; s <- chain[keep, , drop = FALSE]; colnames(s) <- nm
  structure(list(samples = s, acceptance = mean(acc[keep]),
                 mean = stats::setNames(colMeans(s), nm),
                 sd = stats::setNames(apply(s, 2, stats::sd), nm),
                 quantiles = apply(s, 2, stats::quantile, c(0.025, 0.5, 0.975), names = FALSE),
                 pnames = nm), class = "pmcmc_result")
}

#' @export
print.pmcmc_result <- function(x, ...) {
  cat(sprintf("Particle MCMC posterior -- %d samples, acceptance %.2f:\n",
              nrow(x$samples), x$acceptance))
  for (j in seq_along(x$pnames))
    cat(sprintf("  %-10s %.3f  [%.3f, %.3f]\n", x$pnames[j], x$mean[j],
                x$quantiles[1, j], x$quantiles[3, j]))
  invisible(x)
}
