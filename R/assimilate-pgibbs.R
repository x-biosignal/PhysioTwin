# Layer 3 -- data assimilation: particle Gibbs (with ancestor sampling).
#
# particleMCMC() marginalises the latent state and moves the parameters with a
# particle-filter likelihood; particle Gibbs instead SAMPLES the whole state
# trajectory and the parameters in turn. Its state move is a CONDITIONAL SMC: a
# particle filter in which one particle is pinned to the current reference
# trajectory, guaranteeing the sampler cannot collapse. Ancestor sampling
# (Lindsten-Jordan-Schon) re-draws the reference particle's ancestry each step,
# which breaks the path degeneracy of plain particle Gibbs and mixes far better.
# Given a sampled trajectory the parameters have a tractable complete-data
# likelihood, so they are moved by a Metropolis step. Reuses assimilate-filters.R.

# log density of a multivariate normal
.ldmvn <- function(x, mu, S) {
  d <- x - mu; k <- length(x)
  -0.5 * (as.numeric(crossprod(d, solve(S, d))) +
            as.numeric(determinant(S, logarithm = TRUE)$modulus) + k * log(2 * pi))
}

# Conditional SMC with ancestor sampling: sample a new state trajectory given the
# model and a reference trajectory `x_ref` (T x L). Returns a T x L trajectory.
.csmc_as <- function(f, h, y, Q, R, x0, P0, N, x_ref) {
  Y <- as.matrix(y); Tn <- nrow(Y); m <- ncol(Y); L <- length(x0)
  Rinv <- solve(R); LQ <- t(chol(as.matrix(Q))); L0 <- t(chol(as.matrix(P0)))
  Paths <- array(0, dim = c(Tn, L, N)); logw <- numeric(N)
  X1 <- matrix(x0, L, N) + L0 %*% matrix(stats::rnorm(L * N), L, N)
  X1[, N] <- x_ref[1, ]                                     # pin particle N to the reference
  Paths[1, , ] <- X1
  for (i in seq_len(N)) { v <- Y[1, ] - h(X1[, i], 1); logw[i] <- -0.5 * as.numeric(crossprod(v, Rinv %*% v)) }
  for (k in 2:Tn) {
    w <- exp(logw - max(logw)); w <- w / sum(w)
    A <- integer(N)
    A[seq_len(N - 1)] <- sample.int(N, N - 1, replace = TRUE, prob = w)     # ancestors of the free particles
    fprev <- vapply(seq_len(N), function(j) f(Paths[k - 1, , j], k), numeric(L))   # f(x_{k-1}) per particle
    logas <- log(pmax(w, 1e-300)) +                                          # ancestor sampling for the reference
      apply(matrix(fprev, L, N), 2, function(fx) { d <- x_ref[k, ] - fx; -0.5 * as.numeric(crossprod(d, solve(Q, d))) })
    A[N] <- sample.int(N, 1, prob = exp(logas - max(logas)))
    Paths[seq_len(k - 1), , ] <- Paths[seq_len(k - 1), , A, drop = FALSE]    # inherit ancestors' histories
    for (i in seq_len(N - 1)) Paths[k, , i] <- f(Paths[k - 1, , i], k) + LQ %*% stats::rnorm(L)
    Paths[k, , N] <- x_ref[k, ]
    for (i in seq_len(N)) { v <- Y[k, ] - h(Paths[k, , i], k); logw[i] <- -0.5 * as.numeric(crossprod(v, Rinv %*% v)) }
  }
  w <- exp(logw - max(logw)); b <- sample.int(N, 1, prob = w / sum(w))
  matrix(Paths[, , b], Tn, L)
}

#' Particle Gibbs with ancestor sampling
#'
#' Samples the joint posterior of a state-space model's latent trajectory and
#' parameters by particle Gibbs: each iteration draws a new state trajectory with a
#' conditional SMC (ancestor-sampling) sweep given the current parameters, then
#' moves the parameters with a Metropolis step using the tractable complete-data
#' likelihood of that trajectory. Ancestor sampling gives good mixing even for long
#' series where plain particle Gibbs degenerates.
#'
#' @param build_model `function(theta)` returning the state-space model: a list
#'   `f`, `h` (`function(x, k)`), `Q`, `R`, `x0`, `P0`.
#' @param y Observation matrix (`n_steps x m`) or vector.
#' @param prior_logd `function(theta)` log prior density (`-Inf` off the support).
#' @param init Initial parameter vector (its names label the output).
#' @param n_particles Particles in the conditional SMC.
#' @param n_iter Iterations.
#' @param proposal_sd Initial per-parameter proposal SD.
#' @param burn_in Fraction discarded (and over which the proposal adapts).
#' @param seed Optional RNG seed.
#' @return a `pgibbs_result`: parameter posterior `samples`, `acceptance`, `mean`,
#'   `sd`, 2.5/50/97.5% `quantiles`, and the posterior-mean state `trajectory`.
#' @references Andrieu C, Doucet A, Holenstein R (2010) J R Stat Soc B 72:269-342;
#'   Lindsten F, Jordan MI, Schon TB (2014) J Mach Learn Res 15:2145-2184.
#' @seealso [particleMCMC()], [metropolis()], [particleFilter()]
#' @export
#' @examples
#' \donttest{
#' set.seed(1); a <- 0.8; x <- numeric(50); for (k in 2:50) x[k] <- a * x[k-1] + rnorm(1, 0, 0.3)
#' y <- matrix(x + rnorm(50, 0, 0.3), 50, 1)
#' bm <- function(th) list(f = function(x, k) th[1] * x, h = function(x, k) x,
#'                         Q = matrix(0.09), R = matrix(0.09), x0 = 0, P0 = matrix(1))
#' fit <- particleGibbs(bm, y, function(th) dunif(th[1], -1, 1, log = TRUE),
#'                      init = c(a = 0.4), n_particles = 60, n_iter = 600, seed = 1)
#' fit$mean
#' }
particleGibbs <- function(build_model, y, prior_logd, init, n_particles = 100,
                          n_iter = 2000, proposal_sd = 0.1, burn_in = 0.3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Y <- as.matrix(y); Tn <- nrow(Y); d <- length(init)
  nm <- names(init); if (is.null(nm)) nm <- paste0("p", seq_len(d))
  L <- length(build_model(init)$x0)
  cd_ll <- function(theta, x) {                            # complete-data log-likelihood
    m <- build_model(theta)
    ll <- .ldmvn(x[1, ], rep_len(m$x0, L), as.matrix(m$P0))
    for (k in 2:Tn) ll <- ll + .ldmvn(x[k, ], m$f(x[k - 1, ], k), as.matrix(m$Q))
    for (k in seq_len(Tn)) ll <- ll + .ldmvn(Y[k, ], m$h(x[k, ], k), as.matrix(m$R))
    ll
  }
  theta <- as.numeric(init); x_ref <- matrix(0, Tn, L)     # trajectory initialised at 0; PG corrects it
  m0 <- build_model(theta)
  x_ref <- .csmc_as(m0$f, m0$h, Y, m0$Q, m0$R, m0$x0, m0$P0, n_particles, x_ref)
  cdll <- cd_ll(theta, x_ref); lpri <- prior_logd(theta)
  chain <- matrix(0, n_iter, d); acc <- logical(n_iter); traj_sum <- matrix(0, Tn, L)
  Lp <- diag(proposal_sd, d); bwin <- floor(burn_in * n_iter)
  for (i in seq_len(n_iter)) {
    m <- build_model(theta)                                # (1) sample a trajectory given theta
    x_ref <- .csmc_as(m$f, m$h, Y, m$Q, m$R, m$x0, m$P0, n_particles, x_ref)
    cdll <- cd_ll(theta, x_ref)
    thp <- theta + as.numeric(Lp %*% stats::rnorm(d))      # (2) move parameters given the trajectory
    lprip <- prior_logd(thp)
    if (is.finite(lprip)) {
      cdllp <- cd_ll(thp, x_ref)
      if (is.finite(cdllp) && log(stats::runif(1)) < (cdllp + lprip - cdll - lpri)) {
        theta <- thp; cdll <- cdllp; lpri <- lprip; acc[i] <- TRUE
      }
    }
    chain[i, ] <- theta
    if (i > bwin) traj_sum <- traj_sum + x_ref
    if (i >= 100 && i <= bwin && i %% 50 == 0) {
      C <- stats::cov(chain[max(1, i - 400):i, , drop = FALSE]) * (2.4^2 / d) + diag(1e-8, d)
      Lp <- tryCatch(t(chol(C)), error = function(e) Lp)
    }
  }
  keep <- (bwin + 1):n_iter; s <- chain[keep, , drop = FALSE]; colnames(s) <- nm
  structure(list(samples = s, acceptance = mean(acc[keep]),
                 mean = stats::setNames(colMeans(s), nm),
                 sd = stats::setNames(apply(s, 2, stats::sd), nm),
                 quantiles = apply(s, 2, stats::quantile, c(0.025, 0.5, 0.975), names = FALSE),
                 trajectory = traj_sum / length(keep), pnames = nm), class = "pgibbs_result")
}

#' @export
print.pgibbs_result <- function(x, ...) {
  cat(sprintf("Particle Gibbs posterior -- %d samples, acceptance %.2f:\n", nrow(x$samples), x$acceptance))
  for (j in seq_along(x$pnames))
    cat(sprintf("  %-10s %.3f  [%.3f, %.3f]\n", x$pnames[j], x$mean[j],
                x$quantiles[1, j], x$quantiles[3, j]))
  invisible(x)
}
