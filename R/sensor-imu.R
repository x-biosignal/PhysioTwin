# Layer 2 -- the forward IMU sensor model.
#
# A twin only matches reality if the SENSOR is modelled, not just the motion. This
# maps the true limb motion (angle, angular velocity, angular acceleration) to
# what an inertial measurement unit on the segment would actually output:
#   * gyroscope   = angular velocity + bias + white noise + random-walk drift + scale
#   * accelerometer (2-axis, body frame) = specific force (kinematic acceleration of
#     the sensor + gravity), rotated into the segment frame, + bias + noise + scale
# The deterministic part (`.imu_forward`) is also the measurement model used by the
# unscented Kalman filter in Layer 3.

# deterministic forward IMU: state (theta, omega) + angular accel alpha -> (gyro, ax, ay)
.imu_forward <- function(theta, omega, alpha, r, g) {
  # world-frame acceleration of the sensor (tangential + centripetal)
  aw <- r * alpha * c(-sin(theta), cos(theta)) + r * omega^2 * c(-cos(theta), -sin(theta))
  sf <- aw - c(0, -g)                                   # specific force = a - gravity
  R  <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)  # world->body
  c(gyro = omega, ax = as.numeric(R %*% sf)[1], ay = as.numeric(R %*% sf)[2])
}

#' Specify IMU sensor characteristics
#'
#' Bundles the imperfections of an inertial measurement unit -- bias, white noise,
#' random-walk drift and scale error -- so the same sensor can be reused, and
#' randomised across a population (see [generateTrainingData()]).
#'
#' @param gyro_bias,gyro_noise,gyro_drift,gyro_scale Gyroscope constant bias
#'   (rad/s), white-noise SD (rad/s), random-walk drift SD per step (rad/s), and
#'   multiplicative scale error (1 = perfect).
#' @param accel_bias,accel_noise,accel_scale Accelerometer bias (m/s^2, length 1
#'   or 2), white-noise SD (m/s^2) and scale error.
#' @return an `imu_sensor` object.
#' @seealso [imuMeasure()]
#' @export
#' @examples
#' imuSensor(gyro_bias = 0.02, gyro_noise = 0.05)
imuSensor <- function(gyro_bias = 0, gyro_noise = 0.03, gyro_drift = 0,
                      gyro_scale = 1, accel_bias = 0, accel_noise = 0.15,
                      accel_scale = 1) {
  structure(list(gyro_bias = gyro_bias, gyro_noise = gyro_noise,
                 gyro_drift = gyro_drift, gyro_scale = gyro_scale,
                 accel_bias = rep_len(accel_bias, 2), accel_noise = accel_noise,
                 accel_scale = accel_scale), class = "imu_sensor")
}

#' Measure a movement twin with an IMU (forward sensor model)
#'
#' Turns a simulated twin trajectory into a realistic IMU data stream by applying
#' the forward sensor model and the sensor's imperfections.
#'
#' @param trajectory A `twin_trajectory` from [simulateTwin()].
#' @param sensor An `imu_sensor` from [imuSensor()] (default: mild noise).
#' @param seed Optional RNG seed.
#' @param noise If `FALSE`, return the noise-free forward IMU (the ground truth).
#' @return a data frame `time, gyro, ax, ay` (angular velocity rad/s, body-frame
#'   accelerations m/s^2).
#' @seealso [imuSensor()], [personalizeTwin()]
#' @export
#' @examples
#' sim <- simulateTwin(limbTwin(), duration = 2, dt = 0.01)
#' head(imuMeasure(sim, imuSensor(gyro_noise = 0.04), seed = 1))
imuMeasure <- function(trajectory, sensor = imuSensor(), seed = NULL, noise = TRUE) {
  stopifnot(inherits(trajectory, "twin_trajectory"), inherits(sensor, "imu_sensor"))
  p <- trajectory$twin; n <- length(trajectory$time)
  clean <- t(vapply(seq_len(n), function(i)
    .imu_forward(trajectory$theta[i], trajectory$omega[i], trajectory$alpha[i],
                 p$sensor_r, p$g), numeric(3)))
  if (!noise) return(data.frame(time = trajectory$time, gyro = clean[, 1],
                                ax = clean[, 2], ay = clean[, 3]))
  if (!is.null(seed)) set.seed(seed)
  drift <- cumsum(stats::rnorm(n, 0, sensor$gyro_drift))
  gyro <- sensor$gyro_scale * clean[, 1] + sensor$gyro_bias + drift +
    stats::rnorm(n, 0, sensor$gyro_noise)
  ax <- sensor$accel_scale * clean[, 2] + sensor$accel_bias[1] +
    stats::rnorm(n, 0, sensor$accel_noise)
  ay <- sensor$accel_scale * clean[, 3] + sensor$accel_bias[2] +
    stats::rnorm(n, 0, sensor$accel_noise)
  data.frame(time = trajectory$time, gyro = gyro, ax = ax, ay = ay)
}

#' @export
print.imu_sensor <- function(x, ...) {
  cat(sprintf("IMU sensor -- gyro(bias %.3g, noise %.3g, drift %.3g, scale %.3g), accel(noise %.3g)\n",
              x$gyro_bias, x$gyro_noise, x$gyro_drift, x$gyro_scale, x$accel_noise))
  invisible(x)
}
