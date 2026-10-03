# Layer 2 -- forward sensor model: optical motion-capture markers.
#
# A marker-based motion-capture system does not return the true trajectory: each
# marker carries positional noise, and -- whenever the line of sight is blocked --
# drops out for a run of frames (occlusion), leaving gaps that must be filled. This
# maps a clean marker trajectory to what the system records: Gaussian jitter,
# contiguous occlusion gaps (a two-state visible/occluded Markov process, so gaps
# are runs, not isolated frames), optional quantisation and spurious "ghost"
# samples. `fillMarkerGaps()` recovers the gaps by linear interpolation.

#' Specify motion-capture marker characteristics
#'
#' Bundles the imperfections of an optical marker: positional noise, occlusion
#' dropout (contiguous gaps), quantisation and ghost samples.
#'
#' @param noise Gaussian position-noise SD (same units as the trajectory, e.g. m).
#' @param dropout_prob Per-frame probability of STARTING an occlusion while visible.
#' @param dropout_len Mean occlusion length (frames); gaps are geometric.
#' @param ghost_prob Per-frame probability of a spurious (ghost) sample.
#' @param quantization Position resolution (round to this grid; 0 = off).
#' @return a `marker_sensor` object.
#' @seealso [markerMeasure()], [fillMarkerGaps()], [imuSensor()]
#' @export
#' @examples
#' markerSensor(noise = 0.001, dropout_prob = 0.03, dropout_len = 8)
markerSensor <- function(noise = 0.002, dropout_prob = 0.02, dropout_len = 10,
                         ghost_prob = 0, quantization = 0) {
  structure(list(noise = noise, dropout_prob = dropout_prob, dropout_len = dropout_len,
                 ghost_prob = ghost_prob, quantization = quantization),
            class = "marker_sensor")
}

#' Measure a marker trajectory with an optical mocap system
#'
#' Applies the forward marker sensor model to a clean trajectory: adds Gaussian
#' noise, carves contiguous occlusion gaps (visible/occluded Markov process),
#' optionally quantises and injects ghost samples. Occluded frames are returned as
#' `NA` rows.
#'
#' @param position Clean marker trajectory: an `n x d` matrix (or length-`n`
#'   vector) of marker positions over time.
#' @param sensor A `marker_sensor` from [markerSensor()].
#' @param seed Optional RNG seed.
#' @return a `marker_measurement`: `position` (measured, `NA` rows where occluded),
#'   `visible` (logical per frame), `gap_fraction`, and the `sensor`.
#' @seealso [markerSensor()], [fillMarkerGaps()]
#' @export
#' @examples
#' pos <- cbind(sin(seq(0, 6, length.out = 300)), cos(seq(0, 6, length.out = 300)))
#' m <- markerMeasure(pos, markerSensor(dropout_prob = 0.02, dropout_len = 12), seed = 1)
#' m$gap_fraction
markerMeasure <- function(position, sensor = markerSensor(), seed = NULL) {
  stopifnot(inherits(sensor, "marker_sensor"))
  P <- if (is.matrix(position)) position else matrix(position, ncol = 1)
  n <- nrow(P); d <- ncol(P)
  if (!is.null(seed)) set.seed(seed)
  meas <- P + matrix(stats::rnorm(n * d, 0, sensor$noise), n, d)
  if (sensor$quantization > 0) meas <- round(meas / sensor$quantization) * sensor$quantization
  if (sensor$ghost_prob > 0) {                            # spurious large jumps on some frames
    gh <- stats::runif(n) < sensor$ghost_prob
    meas[gh, ] <- meas[gh, ] + matrix(stats::rnorm(sum(gh) * d, 0, 50 * sensor$noise), sum(gh), d)
  }
  visible <- rep(TRUE, n); occ <- FALSE                  # two-state occlusion Markov chain
  for (i in seq_len(n)) {
    visible[i] <- !occ
    if (occ) { if (stats::runif(1) < 1 / sensor$dropout_len) occ <- FALSE }
    else     { if (stats::runif(1) < sensor$dropout_prob)    occ <- TRUE }
  }
  meas[!visible, ] <- NA
  structure(list(position = meas, visible = visible, gap_fraction = mean(!visible),
                 sensor = sensor), class = "marker_measurement")
}

#' Fill occlusion gaps in a measured marker trajectory
#'
#' Linear-interpolates the `NA` runs (occlusions) of a [markerMeasure()] result,
#' column by column, holding the endpoints flat where a gap reaches the boundary.
#'
#' @param m a `marker_measurement` from [markerMeasure()], or a numeric matrix
#'   with `NA` rows.
#' @return the trajectory matrix with gaps interpolated.
#' @seealso [markerMeasure()]
#' @export
#' @examples
#' pos <- matrix(cumsum(rnorm(200)), 200, 1)
#' filled <- fillMarkerGaps(markerMeasure(pos, markerSensor(dropout_prob = 0.05), seed = 1))
fillMarkerGaps <- function(m) {
  P <- if (inherits(m, "marker_measurement")) m$position else as.matrix(m)
  n <- nrow(P); idx <- seq_len(n)
  for (j in seq_len(ncol(P))) {
    ok <- !is.na(P[, j])
    if (sum(ok) >= 2) P[, j] <- stats::approx(idx[ok], P[ok, j], xout = idx, rule = 2)$y
  }
  P
}

#' @export
print.marker_sensor <- function(x, ...) {
  cat(sprintf("Marker sensor -- noise %.3g, dropout(p %.3g, mean len %.3g), ghost %.3g, quant %.3g\n",
              x$noise, x$dropout_prob, x$dropout_len, x$ghost_prob, x$quantization))
  invisible(x)
}

#' @export
print.marker_measurement <- function(x, ...) {
  cat(sprintf("Marker measurement -- %d frames, %.1f%% occluded\n",
              length(x$visible), 100 * x$gap_fraction)); invisible(x)
}
