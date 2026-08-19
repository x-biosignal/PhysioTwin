# Layer 2 -- forward sensor model: a multi-channel surface-EMG array with crosstalk
# and a realistic amplifier chain.
#
# A surface-EMG array does not record each muscle in isolation: neighbouring
# electrodes pick up a fraction of each other's signal (crosstalk), and the
# amplifier chain adds its own character -- mains interference rejected by a finite
# common-mode-rejection ratio (CMRR), a DC-blocking high-pass, a passband and
# thermal noise. This maps per-muscle activation envelopes to the multi-channel
# voltages the array actually records. Band-pass reuses the FFT filter of sensor-emg.R.

#' Electrode crosstalk mixing matrix
#'
#' A symmetric, energy-normalised mixing matrix: each channel keeps most of its own
#' signal and picks up `crosstalk` of an immediate neighbour, `crosstalk^2` of the
#' next, and so on (rows sum to one). `crosstalk = 0` gives the identity.
#'
#' @param channels Number of channels.
#' @param crosstalk Per-neighbour crosstalk fraction (0..1).
#' @return a `channels x channels` mixing matrix.
#' @seealso [emgArray()], [emgArrayMeasure()]
#' @export
#' @examples
#' round(crosstalkMatrix(3, 0.2), 2)
crosstalkMatrix <- function(channels, crosstalk = 0.15) {
  d <- abs(outer(seq_len(channels), seq_len(channels), "-"))
  X <- crosstalk ^ d                       # 1 on the diagonal, crosstalk^distance off it
  X / rowSums(X)                           # energy-preserving: each channel's weights sum to 1
}

# DC-blocking one-pole high-pass: y[n] = a (y[n-1] + x[n] - x[n-1])
.one_pole_hp <- function(x, a) {
  y <- numeric(length(x)); px <- 0; py <- 0
  for (i in seq_along(x)) { y[i] <- a * (py + x[i] - px); px <- x[i]; py <- y[i] }
  y
}

#' Specify a surface-EMG array (crosstalk + amplifier)
#'
#' Bundles a multi-channel surface-EMG measurement chain: inter-electrode crosstalk,
#' amplifier gain, a common-mode-rejected mains interference, a DC-blocking
#' high-pass, a passband and thermal noise.
#'
#' @param channels Number of electrode channels.
#' @param crosstalk Per-neighbour crosstalk fraction (see [crosstalkMatrix()]).
#' @param gain Amplifier gain.
#' @param noise Thermal-noise SD (referred to the electrode).
#' @param powerline_hz,powerline_amp Mains frequency (Hz) and its (common-mode) amplitude.
#' @param cmrr_db Common-mode rejection ratio (dB) attenuating the mains interference.
#' @param hp_cutoff DC-blocking high-pass corner (Hz).
#' @param offset Electrode DC offset (a constant the DC-blocking high-pass removes).
#' @param bandwidth Passband `c(low, high)` (Hz).
#' @param fs Sampling rate (Hz).
#' @return an `emg_array` object.
#' @seealso [emgArrayMeasure()], [crosstalkMatrix()], [emgElectrode()]
#' @export
#' @examples
#' emgArray(channels = 4, crosstalk = 0.2, cmrr_db = 90)
emgArray <- function(channels = 4, crosstalk = 0.15, gain = 1000, noise = 5e-3,
                     powerline_hz = 50, powerline_amp = 0.05, cmrr_db = 80,
                     hp_cutoff = 10, offset = 0, bandwidth = c(20, 450), fs = 1000) {
  structure(list(channels = channels, crosstalk = crosstalk, gain = gain, noise = noise,
                 powerline_hz = powerline_hz, powerline_amp = powerline_amp,
                 cmrr_db = cmrr_db, hp_cutoff = hp_cutoff, offset = offset,
                 bandwidth = bandwidth, fs = fs), class = "emg_array")
}

#' Measure per-muscle activations with a surface-EMG array (forward model)
#'
#' Turns a matrix of per-muscle activation envelopes into the multi-channel raw
#' sEMG an electrode array records: each channel's activation-modulated Gaussian
#' sEMG is band-limited, the channels are MIXED through the crosstalk matrix, and
#' the amplifier adds CMRR-rejected mains interference, a DC-blocking high-pass,
#' thermal noise and gain.
#'
#' @param activations An `n_time x channels` matrix of non-negative activation
#'   envelopes (one column per muscle/channel).
#' @param sensor An `emg_array` from [emgArray()].
#' @param seed Optional RNG seed.
#' @return an `emg_array_measurement`: `raw` (`n_time x channels`), `clean`
#'   (uncontaminated, un-mixed, amplified), `time`, `fs`, the `crosstalk` matrix, sensor.
#' @seealso [emgArray()], [crosstalkMatrix()]
#' @export
#' @examples
#' act <- matrix(pmax(0, sin(outer(seq(0, 6, len = 1200), 1:3))), ncol = 3)
#' m <- emgArrayMeasure(act, emgArray(channels = 3, crosstalk = 0.2), seed = 1)
#' dim(m$raw)
emgArrayMeasure <- function(activations, sensor = emgArray(), seed = NULL) {
  stopifnot(inherits(sensor, "emg_array"))
  A <- if (is.matrix(activations)) activations else matrix(activations, ncol = 1)
  n <- nrow(A); ch <- ncol(A); fs <- sensor$fs; t <- (seq_len(n) - 1) / fs
  if (!is.null(seed)) set.seed(seed)
  clean <- vapply(seq_len(ch), function(j)                       # activation-modulated Gaussian sEMG
    .fft_band(stats::rnorm(n) * A[, j], fs, sensor$bandwidth[1], sensor$bandwidth[2]), numeric(n))
  X <- crosstalkMatrix(ch, sensor$crosstalk)
  mixed <- clean %*% t(X)                                        # measured ch = weighted sum of clean
  pl_amp <- sensor$powerline_amp * 10^(-sensor$cmrr_db / 20)     # common-mode mains after CMRR
  pl <- pl_amp * sin(2 * pi * sensor$powerline_hz * t)
  a <- exp(-2 * pi * sensor$hp_cutoff / fs)                      # DC-blocking pole
  raw <- matrix(0, n, ch)
  for (j in seq_len(ch))                                         # electrode offset is removed by the HP
    raw[, j] <- sensor$gain * .one_pole_hp(mixed[, j] + pl + sensor$offset +
                                             stats::rnorm(n, 0, sensor$noise), a)
  structure(list(raw = raw, clean = sensor$gain * clean, time = t, fs = fs,
                 crosstalk = X, sensor = sensor), class = "emg_array_measurement")
}

#' @export
print.emg_array <- function(x, ...) {
  cat(sprintf("EMG array -- %d channels, crosstalk %.2f, CMRR %g dB, band %g-%g Hz @ %g Hz\n",
              x$channels, x$crosstalk, x$cmrr_db, x$bandwidth[1], x$bandwidth[2], x$fs))
  invisible(x)
}

#' @export
print.emg_array_measurement <- function(x, ...) {
  cat(sprintf("EMG array measurement -- %d samples x %d channels @ %g Hz\n",
              nrow(x$raw), ncol(x$raw), x$fs)); invisible(x)
}
