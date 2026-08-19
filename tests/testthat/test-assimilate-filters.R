# Alternative nonlinear filters: EKF / particle / ensemble Kalman.
# Ground truth is a linear-Gaussian model, where the exact Kalman filter is
# available in closed form -- EKF must match it, PF/EnKF within Monte-Carlo error.

# constant-velocity linear-Gaussian model + its exact Kalman filter
lin_model <- function(seed = 1, n = 60) {
  set.seed(seed)
  dt <- 0.1; A <- matrix(c(1, 0, dt, 1), 2, 2); H <- matrix(c(1, 0), 1, 2)
  Q <- diag(c(1e-3, 1e-2)); R <- matrix(0.05); x0 <- c(0, 1); P0 <- diag(c(1, 1))
  xt <- x0; xs <- matrix(0, n, 2); y <- numeric(n)
  for (k in seq_len(n)) {
    xt <- as.numeric(A %*% xt) + c(rnorm(1, 0, sqrt(Q[1, 1])), rnorm(1, 0, sqrt(Q[2, 2])))
    xs[k, ] <- xt; y[k] <- xt[1] + rnorm(1, 0, sqrt(R))
  }
  list(A = A, H = H, Q = Q, R = R, x0 = x0, P0 = P0, y = matrix(y, n, 1), truth = xs,
       f = function(x, k) as.numeric(A %*% x), h = function(x, k) as.numeric(H %*% x))
}
exact_kf <- function(mo) {
  n <- nrow(mo$y); x <- mo$x0; P <- mo$P0; X <- matrix(0, n, 2)
  for (k in seq_len(n)) {
    xp <- mo$A %*% x; Pp <- mo$A %*% P %*% t(mo$A) + mo$Q
    S <- mo$H %*% Pp %*% t(mo$H) + mo$R; K <- Pp %*% t(mo$H) %*% solve(S)
    x <- xp + K %*% (mo$y[k, ] - mo$H %*% xp); P <- (diag(2) - K %*% mo$H) %*% Pp
    X[k, ] <- x
  }
  X
}

test_that("the extended Kalman filter matches the exact Kalman filter (linear model)", {
  mo <- lin_model(); ref <- exact_kf(mo)
  ekf <- extendedKalmanFilter(mo$f, mo$h, mo$y, mo$Q, mo$R, mo$x0, mo$P0)
  expect_s3_class(ekf, "ekf_result")
  expect_lt(max(abs(ekf$state - ref)), 1e-6)             # EKF is exact when f, h are linear
  expect_output(print(ekf), "EKF")
})

test_that("particle and ensemble filters agree with the Kalman filter within MC error", {
  mo <- lin_model(); ref <- exact_kf(mo)
  pf <- particleFilter(mo$f, mo$h, mo$y, mo$Q, mo$R, mo$x0, mo$P0, n_particles = 4000, seed = 1)
  en <- ensembleKalmanFilter(mo$f, mo$h, mo$y, mo$Q, mo$R, mo$x0, mo$P0, n_ensemble = 300, seed = 1)
  expect_s3_class(pf, "pf_result"); expect_s3_class(en, "enkf_result")
  expect_lt(sqrt(mean((pf$state - ref)^2)), 0.05)
  expect_lt(sqrt(mean((en$state - ref)^2)), 0.06)
  expect_true(all(pf$ess > 0 & pf$ess <= 4000))
})

test_that("all filters track the latent state near the observation-noise floor", {
  mo <- lin_model()
  for (fit in list(extendedKalmanFilter(mo$f, mo$h, mo$y, mo$Q, mo$R, mo$x0, mo$P0),
                   particleFilter(mo$f, mo$h, mo$y, mo$Q, mo$R, mo$x0, mo$P0, n_particles = 2000, seed = 2),
                   ensembleKalmanFilter(mo$f, mo$h, mo$y, mo$Q, mo$R, mo$x0, mo$P0, n_ensemble = 200, seed = 2))) {
    expect_lt(sqrt(mean((fit$state[, 1] - mo$truth[, 1])^2)), 0.3)   # ~ sqrt(R) = 0.22
  }
})

test_that("the EKF accepts an analytic Jacobian", {
  mo <- lin_model()
  ekf_num <- extendedKalmanFilter(mo$f, mo$h, mo$y, mo$Q, mo$R, mo$x0, mo$P0)
  ekf_ana <- extendedKalmanFilter(mo$f, mo$h, mo$y, mo$Q, mo$R, mo$x0, mo$P0,
                                  Fjac = function(x, k) mo$A, Hjac = function(x, k) mo$H)
  expect_lt(max(abs(ekf_num$state - ekf_ana$state)), 1e-4)
})
