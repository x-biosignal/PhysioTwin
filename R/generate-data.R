# Layer 4b -- labelled synthetic training-data generation.
#
# The twin is a labelled data factory: sample physical parameters from a
# population distribution, randomise the sensor (domain randomisation), simulate,
# and emit the sensor stream with the generating parameters as ground-truth
# labels. This produces physically-consistent, richly-labelled training sets for
# machine learning -- the physiological analogue of materials-informatics
# synthetic data -- where every label is a real mechanistic quantity, not a guess.

# a few interpretable features of an IMU stream (ML inputs)
.twin_imu_features <- function(imu, dt) {
  g <- imu$gyro
  sp <- stats::spectrum(g - mean(g), plot = FALSE)
  dom <- sp$freq[which.max(sp$spec)] / dt
  c(gyro_rms = sqrt(mean(g^2)), gyro_range = diff(range(g)),
    ax_rms = sqrt(mean(imu$ax^2)), ay_rms = sqrt(mean(imu$ay^2)),
    gyro_dom_freq = dom)
}

#' Generate a labelled synthetic training set from the twin
#'
#' Samples `n` virtual subjects: physical parameters are drawn from the given
#' ranges, the IMU sensor is randomised (domain randomisation for robustness),
#' each twin is simulated and measured, and the sensor features are returned with
#' the generating parameters as labels. Use it to build ML training data where
#' the labels are true mechanistic quantities.
#'
#' @param n Number of virtual subjects.
#' @param strength,damping,stiffness Length-2 ranges (uniform) for each parameter.
#' @param duration,dt Simulation settings.
#' @param randomize_sensor If `TRUE`, draw a different IMU noise/bias/scale per
#'   subject (domain randomisation).
#' @param features If `TRUE` (default) return summary features; if `FALSE` return
#'   the raw IMU streams alongside the labels.
#' @param seed Optional RNG seed.
#' @return a `twin_dataset`: a data frame of features + label columns
#'   (`strength`, `damping`, `stiffness`), or (when `features = FALSE`) a list
#'   with `labels` and `streams`.
#' @seealso [limbTwin()], [imuMeasure()]
#' @export
#' @examples
#' d <- generateTrainingData(20, duration = 4, dt = 0.02, seed = 1)
#' head(d)
generateTrainingData <- function(n = 100, strength = c(3, 10), damping = c(0.1, 0.8),
                                 stiffness = c(1.5, 4), duration = 5, dt = 0.02,
                                 randomize_sensor = TRUE, features = TRUE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  runif_r <- function(r) stats::runif(1, r[1], r[2])
  rows <- vector("list", n); streams <- vector("list", n)
  for (i in seq_len(n)) {
    p <- list(strength = runif_r(strength), damping = runif_r(damping),
              stiffness = runif_r(stiffness))
    tw <- limbTwin(strength = p$strength, damping = p$damping, stiffness = p$stiffness)
    sim <- simulateTwin(tw, duration, dt, theta0 = 0.1)
    sensor <- if (randomize_sensor)
      imuSensor(gyro_bias = stats::runif(1, -0.03, 0.03),
                gyro_noise = stats::runif(1, 0.01, 0.06),
                gyro_drift = stats::runif(1, 0, 5e-4),
                accel_noise = stats::runif(1, 0.05, 0.25)) else imuSensor()
    imu <- imuMeasure(sim, sensor, noise = TRUE)
    streams[[i]] <- imu
    rows[[i]] <- c(.twin_imu_features(imu, dt), unlist(p))
  }
  labels <- data.frame(do.call(rbind, rows))
  if (features) return(structure(labels, class = c("twin_dataset", "data.frame")))
  structure(list(labels = labels[, c("strength", "damping", "stiffness")],
                 streams = streams), class = "twin_dataset")
}

#' @export
print.twin_dataset <- function(x, ...) {
  if (is.data.frame(x))
    cat(sprintf("Twin dataset -- %d subjects x %d columns (features + labels)\n",
                nrow(x), ncol(x)))
  else cat(sprintf("Twin dataset -- %d subjects (raw streams + labels)\n",
                   nrow(x$labels)))
  invisible(x)
}
