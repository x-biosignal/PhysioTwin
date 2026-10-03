# Layer 2 -- forward sensor model: surface-EMG electrode and amplifier.
#
# A muscle's neural drive is not observed directly; a surface electrode records a
# noisy, band-limited, mains-contaminated voltage. The standard generative model
# of the physiological sEMG is an activation-MODULATED zero-mean Gaussian process
# (its amplitude tracks the activation), which the electrode chain then band-passes
# and contaminates with powerline interference, low-frequency baseline/motion
# wander and thermal noise before amplifying. `emgEnvelope()` recovers the drive
# by rectification and low-pass smoothing. Band-pass and smoothing use base-R FFTs.

# zero out FFT bins outside [lo, hi] Hz (folded to one side) -- an ideal band-pass
.fft_band <- function(x, fs, lo, hi) {
  n <- length(x); f <- (0:(n - 1)) * fs / n; ff <- pmin(f, fs - f)   # fold to [0, fs/2]
  X <- stats::fft(x); X[ff < lo | ff > hi] <- 0
  Re(stats::fft(X, inverse = TRUE)) / n
}

#' Specify surface-EMG electrode/amplifier characteristics
#'
#' Bundles the surface-EMG measurement chain: passband, powerline interference,
#' baseline/motion wander, thermal noise, amplifier gain and clipping.
#'
#' @param gain Amplifier gain.
#' @param noise Thermal (white) noise SD, referred to the electrode.
#' @param powerline_hz Mains frequency (Hz); odd harmonics are added.
#' @param powerline_amp Powerline interference amplitude.
#' @param baseline_amp,baseline_hz Baseline/motion-wander amplitude and frequency (Hz).
#' @param bandwidth Passband `c(low, high)` (Hz).
#' @param fs Sampling rate (Hz).
#' @param clip Saturation level (|output| clipped to this; `Inf` = none).
#' @return an `emg_electrode` object.
#' @seealso [emgMeasure()], [emgEnvelope()], [imuSensor()]
#' @export
#' @examples
#' emgElectrode(powerline_hz = 60, powerline_amp = 0.08)
emgElectrode <- function(gain = 1000, noise = 5e-3, powerline_hz = 50,
                         powerline_amp = 0.05, baseline_amp = 0.02, baseline_hz = 1.0,
                         bandwidth = c(20, 450), fs = 1000, clip = Inf) {
  structure(list(gain = gain, noise = noise, powerline_hz = powerline_hz,
                 powerline_amp = powerline_amp, baseline_amp = baseline_amp,
                 baseline_hz = baseline_hz, bandwidth = bandwidth, fs = fs, clip = clip),
            class = "emg_electrode")
}

#' Measure muscle activation with a surface-EMG electrode (forward model)
#'
#' Turns an activation envelope into a realistic raw sEMG: an activation-modulated
#' Gaussian process band-passed to the electrode's passband (the clean EMG), then
#' contaminated with powerline interference, baseline/motion wander and thermal
#' noise and amplified.
#'
#' @param activation Non-negative activation envelope (0..1), length `n`.
#' @param sensor An `emg_electrode` from [emgElectrode()].
#' @param seed Optional RNG seed.
#' @return an `emg_measurement`: `raw` (contaminated, amplified sEMG), `clean`
#'   (uncontaminated band-limited EMG, amplified), `time`, `activation`, `fs`, `sensor`.
#' @seealso [emgElectrode()], [emgEnvelope()]
#' @export
#' @examples
#' act <- pmax(0, sin(seq(0, 4 * pi, length.out = 2000)))
#' emg <- emgMeasure(act, emgElectrode(fs = 1000), seed = 1)
#' cor(emgEnvelope(emg$clean, 1000), act)
emgMeasure <- function(activation, sensor = emgElectrode(), seed = NULL) {
  stopifnot(inherits(sensor, "emg_electrode"))
  n <- length(activation); fs <- sensor$fs; t <- (seq_len(n) - 1) / fs
  if (!is.null(seed)) set.seed(seed)
  emg <- stats::rnorm(n) * activation                    # activation-modulated Gaussian process
  clean <- .fft_band(emg, fs, sensor$bandwidth[1], sensor$bandwidth[2])
  pl <- sensor$powerline_amp * (sin(2 * pi * sensor$powerline_hz * t) +
        0.3 * sin(2 * pi * 3 * sensor$powerline_hz * t) +
        0.15 * sin(2 * pi * 5 * sensor$powerline_hz * t))          # mains + odd harmonics
  base <- sensor$baseline_amp * sin(2 * pi * sensor$baseline_hz * t)
  raw <- sensor$gain * (clean + pl + base + stats::rnorm(n, 0, sensor$noise))
  raw <- pmax(pmin(raw, sensor$clip), -sensor$clip)
  structure(list(raw = raw, clean = sensor$gain * clean, time = t,
                 activation = activation, fs = fs, sensor = sensor),
            class = "emg_measurement")
}

#' EMG linear envelope (rectify + low-pass)
#'
#' The standard sEMG envelope: full-wave rectification followed by a low-pass at
#' `cutoff` Hz (FFT-based), recovering the activation.
#'
#' @param x A raw or clean EMG signal.
#' @param fs Sampling rate (Hz).
#' @param cutoff Low-pass cut-off (Hz).
#' @return the linear envelope (same length as `x`).
#' @seealso [emgMeasure()]
#' @export
#' @examples
#' emgEnvelope(rnorm(1000) * abs(sin(seq(0, 6, length.out = 1000))), 1000)
emgEnvelope <- function(x, fs, cutoff = 5) {
  n <- length(x); f <- (0:(n - 1)) * fs / n; ff <- pmin(f, fs - f)
  X <- stats::fft(abs(x)); X[ff > cutoff] <- 0            # rectify then low-pass
  Re(stats::fft(X, inverse = TRUE)) / n
}

#' @export
print.emg_electrode <- function(x, ...) {
  cat(sprintf("EMG electrode -- gain %g, band %g-%g Hz, powerline %g Hz (amp %.3g), noise %.3g @ %g Hz\n",
              x$gain, x$bandwidth[1], x$bandwidth[2], x$powerline_hz, x$powerline_amp,
              x$noise, x$fs)); invisible(x)
}

#' @export
print.emg_measurement <- function(x, ...) {
  cat(sprintf("EMG measurement -- %d samples @ %g Hz\n", length(x$raw), x$fs)); invisible(x)
}
