# Layer 1 -- the mechanistic movement twin.
#
# A digital twin starts from a MECHANISTIC model whose parameters mean something
# physical, so that changing them (an intervention) has a principled effect. The
# movement slice uses a parameterised single-joint limb segment: a rigid segment
# rotating about a joint under gravity, passive impedance (stiffness + damping)
# and an active rhythmic drive (the neural command scaled by muscle strength).
# The equation of motion is integrated (RK4) to give the true state trajectory
# that the sensor layer then observes.
#
#   I * theta'' = A * u(t)  -  b * theta'  -  k * (theta - rest)  -  m g l cos(theta)
#
# with parameters that map to clinically meaningful quantities:
#   inertia   I   -- segment morphology
#   strength  A   -- active muscle drive amplitude (what strengthening changes)
#   damping   b   -- passive viscosity
#   stiffness k   -- joint impedance (what spasticity / bracing changes)
#   freq      f   -- movement rhythm; rest -- neutral angle

#' Build a parameterised movement twin (single-joint limb)
#'
#' Constructs a mechanistic single-joint limb model whose parameters map to
#' physical/clinical quantities (segment inertia, active muscle strength, passive
#' damping, joint stiffness). It is the object simulated, personalised and
#' intervened on by the rest of `PhysioTwin`.
#'
#' @param inertia Segment moment of inertia about the joint (kg m^2).
#' @param strength Active-drive amplitude (N m) -- the muscle strength.
#' @param damping Passive viscous damping (N m s / rad).
#' @param stiffness Joint stiffness toward `rest` (N m / rad).
#' @param freq Rhythm of the active drive (Hz).
#' @param rest Neutral joint angle (rad).
#' @param mass,com,g Segment mass (kg), centre-of-mass distance from the joint
#'   (m) and gravity (m/s^2); set `g = 0` for a horizontal-plane movement.
#' @param sensor_r Distance of the IMU from the joint along the segment (m), used
#'   by the sensor layer.
#' @return a `movement_twin` object (a named parameter list with class).
#' @seealso [simulateTwin()], [personalizeTwin()], [insilicoIntervention()]
#' @export
#' @examples
#' tw <- limbTwin(strength = 6, damping = 0.4, stiffness = 3)
#' sim <- simulateTwin(tw, duration = 4, dt = 0.005)
#' range(sim$theta)
limbTwin <- function(inertia = 0.15, strength = 5, damping = 0.3, stiffness = 2.5,
                     freq = 1.0, rest = 0, mass = 2.5, com = 0.2, g = 9.81,
                     sensor_r = 0.25) {
  p <- list(inertia = inertia, strength = strength, damping = damping,
            stiffness = stiffness, freq = freq, rest = rest, mass = mass,
            com = com, g = g, sensor_r = sensor_r)
  stopifnot(p$inertia > 0)
  structure(p, class = "movement_twin")
}

# active drive command u(t) in [-1, 1]
.twin_drive <- function(t, p) sin(2 * pi * p$freq * t)

# state derivative d[theta, omega]/dt
.twin_deriv <- function(t, s, p) {
  theta <- s[1]; omega <- s[2]
  tau <- p$strength * .twin_drive(t, p) - p$damping * omega -
    p$stiffness * (theta - p$rest) - p$mass * p$g * p$com * cos(theta)
  c(omega, tau / p$inertia)
}

#' Simulate a movement twin (forward dynamics)
#'
#' Integrates the twin's equation of motion with a fixed-step 4th-order
#' Runge-Kutta scheme, returning the true joint state trajectory (angle, angular
#' velocity and acceleration).
#'
#' @param twin A `movement_twin` from [limbTwin()].
#' @param duration Simulation time (s).
#' @param dt Time step (s).
#' @param theta0,omega0 Initial angle and angular velocity.
#' @return a `twin_trajectory` list: `time`, `theta`, `omega`, `alpha` (angular
#'   acceleration), and `twin`.
#' @seealso [limbTwin()], [imuMeasure()]
#' @export
#' @examples
#' simulateTwin(limbTwin(), duration = 2, dt = 0.005)$theta[1:5]
simulateTwin <- function(twin, duration = 5, dt = 0.005, theta0 = 0, omega0 = 0) {
  stopifnot(inherits(twin, "movement_twin"))
  n <- as.integer(round(duration / dt)) + 1L
  th <- om <- al <- numeric(n); tt <- (seq_len(n) - 1) * dt
  s <- c(theta0, omega0)
  th[1] <- s[1]; om[1] <- s[2]; al[1] <- .twin_deriv(0, s, twin)[2]
  for (i in 2:n) {
    t <- tt[i - 1]
    k1 <- .twin_deriv(t, s, twin)
    k2 <- .twin_deriv(t + dt / 2, s + dt / 2 * k1, twin)
    k3 <- .twin_deriv(t + dt / 2, s + dt / 2 * k2, twin)
    k4 <- .twin_deriv(t + dt, s + dt * k3, twin)
    s <- s + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    th[i] <- s[1]; om[i] <- s[2]; al[i] <- .twin_deriv(tt[i], s, twin)[2]
  }
  structure(list(time = tt, theta = th, omega = om, alpha = al, twin = twin, dt = dt),
            class = "twin_trajectory")
}

#' @export
print.movement_twin <- function(x, ...) {
  cat("Movement twin (single-joint limb)\n")
  cat(sprintf("  inertia %.3f  strength %.2f  damping %.2f  stiffness %.2f  freq %.2f Hz\n",
              x$inertia, x$strength, x$damping, x$stiffness, x$freq))
  invisible(x)
}

#' @export
print.twin_trajectory <- function(x, ...) {
  cat(sprintf("Twin trajectory -- %.2f s @ dt %g (%d samples)\n",
              max(x$time), x$dt, length(x$time)))
  cat(sprintf("  angle range [%.2f, %.2f] rad, peak |omega| %.2f rad/s\n",
              min(x$theta), max(x$theta), max(abs(x$omega))))
  invisible(x)
}
