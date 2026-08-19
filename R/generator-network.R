# Coordination / network generators: coupled oscillators and a central pattern
# generator.
#
# Physiological rhythms are produced by NETWORKS of oscillators, not single
# templates. Two canonical mechanistic generators:
#   * kuramoto()     -- N phase oscillators with all-to-all coupling; below a
#                       critical coupling they are incoherent, above it they
#                       synchronise (the order parameter jumps) -- the archetype of
#                       collective rhythm / inter-limb / neural synchronisation.
#   * cpgMatsuoka()  -- a Matsuoka half-centre oscillator: two mutually-inhibiting
#                       adapting neurons that produce a self-sustained rhythm -- the
#                       building block of central pattern generators for locomotion.
# Dependency-free base R, RK4 integration.

#' Kuramoto coupled-oscillator network
#'
#' Simulates `n` phase oscillators with natural frequencies `omega` and all-to-all
#' coupling `K` (`theta_i' = omega_i + K/n sum_j sin(theta_j - theta_i)`). The
#' order parameter measures collective synchronisation: near 0 (incoherent) below
#' the critical coupling, near 1 (synchronised) above it.
#'
#' @param n Number of oscillators.
#' @param K Coupling strength.
#' @param omega Natural frequencies (length `n`); default `rnorm(n)`.
#' @param duration,dt Simulation length and step (s).
#' @param seed Optional RNG seed.
#' @return a `kuramoto` object: `time`, `phase` (`nt x n`), `order` (order
#'   parameter per step), `signals` (`cos(phase)`), `mean_order`.
#' @references Kuramoto Y (1975); Strogatz SH (2000) Physica D 143:1-20.
#' @seealso [cpgMatsuoka()], [ecgsyn()], [jansenRit()]
#' @export
#' @examples
#' kuramoto(n = 15, K = 3, duration = 20, seed = 1)$mean_order   # synchronised
#' kuramoto(n = 15, K = 0, duration = 20, seed = 1)$mean_order   # incoherent
kuramoto <- function(n = 20, K = 1, omega = NULL, duration = 20, dt = 0.01, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(omega)) omega <- stats::rnorm(n)
  theta <- stats::runif(n, 0, 2 * pi)
  nt <- round(duration / dt)
  Theta <- matrix(0, nt, n); ord <- numeric(nt)
  deriv <- function(th) omega - (K / n) * rowSums(sin(outer(th, th, `-`)))
  for (i in seq_len(nt)) {
    k1 <- deriv(theta); k2 <- deriv(theta + dt/2 * k1)
    k3 <- deriv(theta + dt/2 * k2); k4 <- deriv(theta + dt * k3)
    theta <- theta + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
    Theta[i, ] <- theta; ord[i] <- Mod(mean(exp(1i * theta)))
  }
  keep <- floor(nt / 2):nt                              # after transient
  structure(list(time = (seq_len(nt)) * dt, phase = Theta, order = ord,
                 signals = cos(Theta), mean_order = mean(ord[keep]),
                 n = n, K = K), class = "kuramoto")
}

#' @export
print.kuramoto <- function(x, ...) {
  cat(sprintf("Kuramoto network -- %d oscillators, K = %.2f, order parameter %.2f (%s)\n",
              x$n, x$K, x$mean_order,
              if (x$mean_order > 0.6) "synchronised" else if (x$mean_order < 0.3) "incoherent" else "partial"))
  invisible(x)
}

#' Matsuoka central pattern generator (half-centre oscillator)
#'
#' Simulates a Matsuoka two-neuron half-centre oscillator: two mutually-inhibiting
#' neurons with fatigue/adaptation produce a self-sustained anti-phase rhythm --
#' the mechanistic building block of locomotor central pattern generators. The
#' output is the difference of the two rectified activities.
#'
#' @param tau,tau2 Rise and adaptation time constants (s).
#' @param beta Adaptation (self-inhibition) strength.
#' @param gamma Mutual-inhibition strength.
#' @param u Tonic drive.
#' @param duration,dt Simulation length and step (s).
#' @return a `cpg` object: `time`, `output` (rhythmic drive `y1 - y2`), `y1`, `y2`
#'   (the two neuron outputs), and `frequency` (Hz, estimated from the output).
#' @references Matsuoka K (1985) Biol Cybern 52:367-376.
#' @seealso [kuramoto()], [jansenRit()]
#' @export
#' @examples
#' cpg <- cpgMatsuoka(duration = 20)
#' cpg$frequency        # a sustained locomotor-like rhythm
cpgMatsuoka <- function(tau = 0.25, tau2 = 0.5, beta = 2.5, gamma = 2.5, u = 1,
                        duration = 20, dt = 0.005) {
  relu <- function(x) pmax(x, 0)
  deriv <- function(s) {
    x1 <- s[1]; v1 <- s[2]; x2 <- s[3]; v2 <- s[4]
    c((-x1 - beta * v1 - gamma * relu(x2) + u) / tau,
      (-v1 + relu(x1)) / tau2,
      (-x2 - beta * v2 - gamma * relu(x1) + u) / tau,
      (-v2 + relu(x2)) / tau2)
  }
  nt <- round(duration / dt); s <- c(0.2, 0, -0.2, 0)   # asymmetric init breaks symmetry
  Y1 <- Y2 <- numeric(nt)
  for (i in seq_len(nt)) {
    k1 <- deriv(s); k2 <- deriv(s + dt/2 * k1)
    k3 <- deriv(s + dt/2 * k2); k4 <- deriv(s + dt * k3)
    s <- s + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
    Y1[i] <- relu(s[1]); Y2[i] <- relu(s[3])
  }
  y <- Y1 - Y2; t <- seq_len(nt) * dt
  keep <- floor(nt / 3):nt; yk <- y[keep] - mean(y[keep])
  sp <- stats::spectrum(yk, plot = FALSE); freq <- sp$freq[which.max(sp$spec)] / dt
  structure(list(time = t, output = y, y1 = Y1, y2 = Y2, frequency = freq),
            class = "cpg")
}

#' @export
print.cpg <- function(x, ...) {
  cat(sprintf("Matsuoka CPG -- %.2f s, rhythm %.2f Hz, output range [%.2f, %.2f]\n",
              max(x$time), x$frequency, min(x$output), max(x$output)))
  invisible(x)
}
