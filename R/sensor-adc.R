# Layer 2 -- forward sensor model: analogue-to-digital conversion artefacts.
#
# Digitising an analogue biosignal quantises its amplitude to a finite bit depth
# and samples it at a finite rate, possibly with timing jitter. These helpers add
# those converter artefacts -- amplitude quantisation (with optional dither),
# resampling with sampling jitter, and the theoretical/empirical quantisation
# noise -- so a simulated signal can be degraded the way a real data-logger would.

#' Quantise a signal to a finite bit depth
#'
#' Rounds a signal to `bits` resolution over `range`, as an analogue-to-digital
#' converter does. Optional (subtractive) dither randomises the quantisation error,
#' decorrelating it from the signal.
#'
#' @param signal Numeric signal.
#' @param bits Bit depth (e.g. 8, 12, 16).
#' @param range `c(lo, hi)` full-scale range; defaults to the signal's range.
#' @param dither If `TRUE`, add uniform dither of one LSB peak-to-peak before rounding.
#' @return the quantised signal (clipped to `range`).
#' @seealso [adcSample()], [quantizationNoise()]
#' @export
#' @examples
#' q <- adcQuantize(sin(seq(0, 6, length.out = 500)), bits = 8)
#' length(unique(round(q, 6))) <= 256
adcQuantize <- function(signal, bits = 12, range = NULL, dither = FALSE) {
  if (is.null(range)) range <- range(signal)
  lo <- range[1]; lsb <- diff(range) / (2^bits - 1)
  x <- if (dither) signal + stats::runif(length(signal), -lsb / 2, lsb / 2) else signal
  q <- round((x - lo) / lsb) * lsb + lo
  pmin(pmax(q, range[1]), range[2])
}

#' Resample a signal to a new rate with optional sampling jitter
#'
#' Resamples `signal` from `fs_in` to `fs_out` by linear interpolation, optionally
#' perturbing each output sample time by Gaussian timing `jitter` -- the aperture
#' jitter of a real converter's clock.
#'
#' @param signal Numeric signal sampled at `fs_in`.
#' @param fs_in,fs_out Input and output sampling rates (Hz).
#' @param jitter Sampling-time jitter SD, in input samples.
#' @return a list: `signal` (resampled), `fs` (`= fs_out`), `time`.
#' @seealso [adcQuantize()]
#' @export
#' @examples
#' down <- adcSample(sin(seq(0, 12, length.out = 2000)), fs_in = 1000, fs_out = 250)
#' down$fs
adcSample <- function(signal, fs_in, fs_out, jitter = 0) {
  n <- length(signal); dur <- (n - 1) / fs_in
  t_out <- seq(0, dur, by = 1 / fs_out)
  tt <- if (jitter > 0) t_out + stats::rnorm(length(t_out), 0, jitter / fs_in) else t_out
  y <- stats::approx((seq_len(n) - 1) / fs_in, signal, xout = tt, rule = 2)$y
  list(signal = y, fs = fs_out, time = t_out)
}

#' Quantisation noise of an analogue-to-digital conversion
#'
#' The theoretical quantisation noise (least-significant-bit squared over twelve)
#' and the empirical error of quantising `signal` to `bits` over `range`.
#'
#' @param signal Numeric signal.
#' @param bits Bit depth.
#' @param range `c(lo, hi)` full-scale range; defaults to the signal's range.
#' @return a named numeric vector: `lsb`, `theoretical_sd` (`= lsb/sqrt(12)`),
#'   `empirical_sd`, `theoretical_power` (`= lsb^2/12`).
#' @seealso [adcQuantize()]
#' @export
#' @examples
#' quantizationNoise(sin(seq(0, 20, length.out = 4000)), bits = 12)
quantizationNoise <- function(signal, bits, range = NULL) {
  if (is.null(range)) range <- range(signal)
  lsb <- diff(range) / (2^bits - 1)
  c(lsb = lsb, theoretical_sd = lsb / sqrt(12),
    empirical_sd = stats::sd(signal - adcQuantize(signal, bits, range)),
    theoretical_power = lsb^2 / 12)
}
