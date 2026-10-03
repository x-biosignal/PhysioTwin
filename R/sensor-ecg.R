# Layer 2 -- forward sensor model: ECG lead / electrode.
#
# A recorded ECG lead is the clean cardiac signal plus the electrode chain's
# artefacts: low-frequency BASELINE WANDER (respiration + electrode half-cell
# drift), 50/60 Hz POWERLINE interference, high-frequency muscle (EMG) noise, and
# optional motion artefact, with an electrode gain and offset. This maps a clean
# ECG (e.g. from ecgsyn()) to the contaminated lead a monitor would record.
# `detectRpeaks()` finds the R waves so recovery can be checked.

#' Specify ECG lead / electrode characteristics
#'
#' Bundles the ECG recording chain's artefacts: baseline wander, powerline
#' interference, muscle noise, motion artefact, gain and offset.
#'
#' @param gain Electrode/amplifier gain.
#' @param noise High-frequency muscle-noise SD (mV).
#' @param powerline_hz Mains frequency (Hz).
#' @param powerline_amp Powerline interference amplitude (mV).
#' @param baseline_amp,baseline_hz Baseline-wander amplitude (mV) and frequency (Hz,
#'   typically respiratory ~0.25 Hz).
#' @param motion_amp Motion-artefact amplitude (mV; a slow random walk).
#' @param offset Electrode DC offset (mV).
#' @param fs Sampling rate (Hz).
#' @return an `ecg_lead` object.
#' @seealso [ecgMeasure()], [detectRpeaks()], [ecgsyn()]
#' @export
#' @examples
#' ecgLead(baseline_amp = 0.2, powerline_hz = 60)
ecgLead <- function(gain = 1, noise = 0.01, powerline_hz = 50, powerline_amp = 0.05,
                    baseline_amp = 0.15, baseline_hz = 0.25, motion_amp = 0,
                    offset = 0, fs = 250) {
  structure(list(gain = gain, noise = noise, powerline_hz = powerline_hz,
                 powerline_amp = powerline_amp, baseline_amp = baseline_amp,
                 baseline_hz = baseline_hz, motion_amp = motion_amp, offset = offset,
                 fs = fs), class = "ecg_lead")
}

#' Measure a clean ECG through a lead/electrode (forward model)
#'
#' Adds the electrode chain's artefacts -- baseline wander, powerline interference,
#' muscle noise, optional motion artefact, gain and offset -- to a clean ECG.
#'
#' @param ecg Clean ECG signal (mV), e.g. from [ecgsyn()].
#' @param sensor An `ecg_lead` from [ecgLead()].
#' @param seed Optional RNG seed.
#' @return an `ecg_measurement`: `signal` (recorded lead, mV), `clean`, `time`,
#'   `fs`, `sensor`.
#' @seealso [ecgLead()], [detectRpeaks()]
#' @export
#' @examples
#' e <- ecgsyn(duration = 10, seed = 1)
#' rec <- ecgMeasure(e$ecg, ecgLead(baseline_amp = 0.2, fs = e$sfecg), seed = 1)
#' sd(rec$signal)
ecgMeasure <- function(ecg, sensor = ecgLead(), seed = NULL) {
  stopifnot(inherits(sensor, "ecg_lead"))
  n <- length(ecg); fs <- sensor$fs; t <- (seq_len(n) - 1) / fs
  if (!is.null(seed)) set.seed(seed)
  base <- sensor$baseline_amp * sin(2 * pi * sensor$baseline_hz * t)
  pl <- sensor$powerline_amp * sin(2 * pi * sensor$powerline_hz * t)
  motion <- if (sensor$motion_amp > 0)
    sensor$motion_amp * cumsum(stats::rnorm(n)) / sqrt(n) else 0
  sig <- sensor$gain * ecg + base + pl + motion +
    stats::rnorm(n, 0, sensor$noise) + sensor$offset
  structure(list(signal = sig, clean = ecg, time = t, fs = fs, sensor = sensor),
            class = "ecg_measurement")
}

#' Remove baseline wander and detect R-peaks in an ECG
#'
#' High-passes the signal (subtracting a ~1 s moving-average trend) and returns the
#' indices of R waves -- local maxima above a fraction of the signal range,
#' separated by a physiological refractory gap.
#'
#' @param x An ECG signal (numeric).
#' @param fs Sampling rate (Hz).
#' @param refractory Minimum spacing between R-peaks (s).
#' @param thresh Peak threshold as a fraction of the (de-trended) signal range.
#' @return the integer indices of detected R-peaks.
#' @seealso [ecgMeasure()]
#' @export
#' @examples
#' e <- ecgsyn(duration = 10, seed = 1)
#' length(detectRpeaks(e$ecg, e$sfecg))
detectRpeaks <- function(x, fs, refractory = 0.3, thresh = 0.5) {
  w <- max(3, round(fs))                                  # ~1 s moving-average trend
  trend <- stats::filter(x, rep(1 / w, w), sides = 2)
  hp <- x - ifelse(is.na(trend), stats::median(x, na.rm = TRUE), trend)  # baseline removed
  sw <- max(1L, round(fs / 40))                          # smooth into the QRS band (kills
  if (sw > 1) {                                          # powerline + muscle noise)
    hp <- as.numeric(stats::filter(hp, rep(1 / sw, sw), sides = 2)); hp[is.na(hp)] <- 0
  }
  lvl <- thresh * max(hp)                                # R-peaks are the tall positive deflections
  gap <- round(refractory * fs); n <- length(hp); peaks <- integer(0); i <- 2L
  while (i < n) {
    if (hp[i] > lvl && hp[i] >= hp[i - 1] && hp[i] > hp[i + 1]) {
      peaks <- c(peaks, i); i <- i + gap
    } else i <- i + 1L
  }
  peaks
}

#' @export
print.ecg_lead <- function(x, ...) {
  cat(sprintf("ECG lead -- gain %g, baseline %.3g mV @ %g Hz, powerline %g Hz (amp %.3g), noise %.3g\n",
              x$gain, x$baseline_amp, x$baseline_hz, x$powerline_hz, x$powerline_amp,
              x$noise)); invisible(x)
}

#' @export
print.ecg_measurement <- function(x, ...) {
  cat(sprintf("ECG measurement -- %d samples @ %g Hz\n", length(x$signal), x$fs)); invisible(x)
}
