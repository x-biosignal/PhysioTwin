# Layer 3 -- data assimilation: full Bayesian posteriors.
#
# abcCalibration() gives a rejection-ABC posterior; two stronger engines here:
#   * metropolis()  -- adaptive random-walk Metropolis-Hastings for any tractable
#                      log-posterior (likelihood x prior); the proposal covariance
#                      adapts to the chain during burn-in (Haario), so correlated
#                      posteriors mix well.
#   * abcSMC()       -- sequential Monte-Carlo ABC (Toni et al. 2009): a particle
#                      population is driven through a schedule of shrinking
#                      tolerances with a perturbation kernel and importance weights,
#                      concentrating on the posterior far more efficiently than
#                      one-shot rejection ABC.
# Dependency-free base R.

#' Adaptive random-walk Metropolis-Hastings sampler
#'
#' Samples from a log-posterior by random-walk Metropolis-Hastings. The Gaussian
#' proposal covariance adapts to the empirical covariance of the chain during the
#' burn-in (Haario et al. 2001), then is frozen, so correlated targets are sampled
#' efficiently while the stationary distribution is preserved.
#'
#' @param logpost `function(theta)` returning the log-posterior (up to a constant);
#'   may return `-Inf` outside the support.
#' @param init Initial parameter vector (its names label the output).
#' @param n_iter Number of iterations.
#' @param proposal_sd Initial per-parameter proposal SD.
#' @param burn_in Fraction discarded (and over which the proposal adapts).
#' @param seed Optional RNG seed.
#' @return an `mcmc_result`: posterior `samples`, `acceptance` rate, `mean`, `sd`,
#'   and 2.5/50/97.5% `quantiles`.
#' @references Metropolis et al. (1953); Haario H, Saksman E, Tamminen J (2001)
#'   Bernoulli 7:223-242 (adaptive Metropolis).
#' @seealso [abcSMC()], [abcCalibration()], [gpCalibrate()]
#' @export
#' @examples
#' # posterior of a normal mean: 20 obs, N(0,10) prior
#' set.seed(1); y <- rnorm(20, 2, 1)
#' lp <- function(m) sum(dnorm(y, m, 1, log = TRUE)) + dnorm(m, 0, 10, log = TRUE)
#' fit <- metropolis(lp, init = c(mu = 0), n_iter = 4000, seed = 1)
#' fit$mean
metropolis <- function(logpost, init, n_iter = 6000, proposal_sd = 0.5,
                       burn_in = 0.4, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  d <- length(init); nm <- names(init); if (is.null(nm)) nm <- paste0("p", seq_len(d))
  L <- diag(proposal_sd, d)                       # proposal Cholesky factor (diagonal at first)
  x <- as.numeric(init); lp <- logpost(x)
  chain <- matrix(0, n_iter, d); acc <- logical(n_iter)
  bwin <- floor(burn_in * n_iter); sd_scale <- 2.4^2 / d
  for (i in seq_len(n_iter)) {
    xp <- x + as.numeric(L %*% stats::rnorm(d)); lpp <- logpost(xp)
    if (is.finite(lpp) && log(stats::runif(1)) < (lpp - lp)) { x <- xp; lp <- lpp; acc[i] <- TRUE }
    chain[i, ] <- x
    if (i >= 200 && i <= bwin && i %% 50 == 0) {   # adapt proposal to the chain (burn-in only)
      C <- stats::cov(chain[max(1, i - 800):i, , drop = FALSE]) * sd_scale + diag(1e-8, d)
      L <- tryCatch(t(chol(C)), error = function(e) L)
    }
  }
  keep <- (bwin + 1):n_iter; s <- chain[keep, , drop = FALSE]; colnames(s) <- nm
  structure(list(samples = s, acceptance = mean(acc[keep]),
                 mean = stats::setNames(colMeans(s), nm),
                 sd = stats::setNames(apply(s, 2, stats::sd), nm),
                 quantiles = apply(s, 2, stats::quantile, c(0.025, 0.5, 0.975), names = FALSE),
                 pnames = nm), class = "mcmc_result")
}

#' @export
print.mcmc_result <- function(x, ...) {
  cat(sprintf("MCMC posterior -- %d samples, acceptance %.2f:\n", nrow(x$samples), x$acceptance))
  for (j in seq_along(x$pnames))
    cat(sprintf("  %-10s %.3f  [%.3f, %.3f]\n", x$pnames[j], x$mean[j],
                x$quantiles[1, j], x$quantiles[3, j]))
  invisible(x)
}

# density of x under N(rows of M, Sigma), one value per row of M
.dmvn_rows <- function(x, M, Sigma) {
  p <- length(x); Si <- solve(Sigma)
  logdet <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
  D <- sweep(M, 2, x, "-"); q <- rowSums((D %*% Si) * D)
  as.numeric(exp(-0.5 * (q + logdet + p * log(2 * pi))))
}

#' Sequential Monte-Carlo ABC (likelihood-free posterior by tolerance annealing)
#'
#' A sequential Monte-Carlo approximate Bayesian computation sampler. A population
#' of parameter particles is carried through a schedule of shrinking tolerances:
#' at each level a particle is drawn from the previous population, perturbed by a
#' Gaussian kernel (covariance twice the population's), simulated, and accepted if
#' its summary distance is below the current tolerance; importance weights correct
#' for the kernel. The tolerance for each level is an adaptive quantile of the
#' previous population's distances. This concentrates on the posterior with far
#' fewer simulations than one-shot rejection ABC.
#'
#' @param simulate `function(theta)` returning summary statistics.
#' @param observed Observed summary statistics.
#' @param prior `function()` returning one prior draw (named).
#' @param prior_density Optional `function(theta)` prior density (for the weights);
#'   `NULL` assumes a (bounded) uniform prior.
#' @param n_particles Particles per population.
#' @param n_populations Number of tolerance levels.
#' @param quantile Tolerance schedule: each level's tolerance is this quantile of
#'   the previous population's distances.
#' @param max_tries Per-particle acceptance attempts before giving up (guards a
#'   too-tight tolerance).
#' @param seed Optional RNG seed.
#' @return an `abc_smc`: weighted posterior `samples`, `weights`, the `tolerance`
#'   schedule, and the weighted `mean`/`sd`.
#' @references Toni T et al. (2009) J R Soc Interface 6:187-202; Beaumont MA et al.
#'   (2009) Biometrika 96:983-990.
#' @seealso [abcCalibration()], [metropolis()]
#' @export
#' @examples
#' set.seed(1); obs <- c(mean = 2, sd = 1.5)
#' sim <- function(th) { x <- rnorm(200, th[1], th[2]); c(mean(x), sd(x)) }
#' pri <- function() c(m = runif(1, -3, 7), s = runif(1, 0.2, 4))
#' post <- abcSMC(sim, obs, pri, n_particles = 100, n_populations = 4, seed = 1)
#' post$mean
abcSMC <- function(simulate, observed, prior, prior_density = NULL,
                   n_particles = 200, n_populations = 5, quantile = 0.5,
                   max_tries = 2000L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  observed <- as.numeric(observed)
  th0 <- prior(); p <- length(th0); pnames <- names(th0)
  if (is.null(pnames)) pnames <- paste0("theta", seq_len(p))
  if (is.null(prior_density)) prior_density <- function(th) 1
  pre <- t(vapply(seq_len(max(60L, n_particles)),
                  function(i) as.numeric(simulate(prior())), numeric(length(observed))))
  sdev <- apply(pre, 2, stats::sd); sdev[sdev < 1e-8] <- 1
  dfun <- function(th) {                                # invalid proposals (out of support) -> rejected
    s <- tryCatch(suppressWarnings(as.numeric(simulate(th))), error = function(e) NA_real_)
    if (length(s) != length(observed) || anyNA(s) || any(!is.finite(s))) return(Inf)
    sqrt(sum(((s - observed) / sdev)^2))
  }
  part <- matrix(0, n_particles, p, dimnames = list(NULL, pnames)); dvec <- numeric(n_particles)
  for (i in seq_len(n_particles)) { th <- prior(); part[i, ] <- th; dvec[i] <- dfun(th) }
  w <- rep(1 / n_particles, n_particles); tol_sched <- numeric(n_populations)
  tol_sched[1] <- stats::quantile(dvec, quantile, names = FALSE)
  for (t in 2:n_populations) {
    tol <- stats::quantile(dvec, quantile, names = FALSE); tol_sched[t] <- tol
    cw <- stats::cov.wt(part, wt = w, method = "ML"); K <- 2 * cw$cov + diag(1e-9, p)
    Lk <- t(chol(K))
    np <- matrix(0, n_particles, p, dimnames = list(NULL, pnames)); nd <- numeric(n_particles)
    nw <- numeric(n_particles)
    for (i in seq_len(n_particles)) {
      thp <- NULL
      for (tr in seq_len(max_tries)) {
        j <- sample.int(n_particles, 1, prob = w)
        cand <- part[j, ] + as.numeric(Lk %*% stats::rnorm(p))
        if (prior_density(cand) <= 0) next
        d <- dfun(cand)
        if (d < tol) { thp <- cand; break }
      }
      if (is.null(thp)) { thp <- part[j, ]; d <- dvec[j] }     # fallback: keep a prior particle
      np[i, ] <- thp; nd[i] <- d
      nw[i] <- prior_density(thp) / sum(w * .dmvn_rows(thp, part, K))
    }
    part <- np; dvec <- nd; w <- nw / sum(nw)
  }
  cw <- stats::cov.wt(part, wt = w, method = "ML")
  structure(list(samples = part, weights = w, tolerance = tol_sched,
                 mean = stats::setNames(cw$center, pnames),
                 sd = stats::setNames(sqrt(diag(cw$cov)), pnames), pnames = pnames),
            class = "abc_smc")
}

#' @export
print.abc_smc <- function(x, ...) {
  cat(sprintf("SMC-ABC posterior -- %d particles, tolerance %.3g -> %.3g:\n",
              nrow(x$samples), x$tolerance[1], x$tolerance[length(x$tolerance)]))
  for (j in seq_along(x$pnames))
    cat(sprintf("  %-10s %.3f (sd %.3f)\n", x$pnames[j], x$mean[j], x$sd[j]))
  invisible(x)
}
