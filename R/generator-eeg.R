# Physiological EEG generator: the Jansen-Rit neural-mass model.
#
# EEG is not a template + noise; it is the emergent output of a cortical column
# of interacting neural populations. The Jansen-Rit model (Jansen & Rit 1995) is
# the canonical mechanistic generator: a single cortical column of THREE
# populations -- pyramidal cells, excitatory interneurons and inhibitory
# interneurons -- each converting an incoming firing rate to an average
# post-synaptic potential (a second-order "alpha" kernel) and back through a
# sigmoid. With standard parameters the column self-organises into an ~10 Hz
# ALPHA rhythm; changing the excitation/inhibition balance or the external drive
# moves it between resting, alpha and faster/spiking regimes. The scalp EEG is
# the pyramidal post-synaptic potential y1 - y2.
#
# Dependency-free base R, RK4 integration (internally oversampled for accuracy).

# Map a dominant frequency (Hz) to the conventional EEG band.
.eeg_band <- function(hz) {
  if (hz < 0.5) "sub-delta" else if (hz < 4) "delta" else if (hz < 8) "theta" else
    if (hz < 13) "alpha" else if (hz < 30) "beta" else "gamma"
}

# Smoothed, detrended one-sided PSD of the post-transient EEG (Hz vs power).
# A raw periodogram is far too noisy to read a dominant rhythm from -- the peak
# gets pulled onto spurious low-frequency spikes -- so smooth and detrend.
.eeg_psd <- function(eeg, sfeeg) {
  n <- length(eeg); keep <- floor(n / 3):n                # drop transient
  xk <- eeg[keep] - mean(eeg[keep])
  sp <- stats::spec.pgram(xk, spans = c(7, 7), taper = 0.1,
                          detrend = TRUE, plot = FALSE)
  list(freq = sp$freq * sfeeg, spec = sp$spec)
}

#' Aperiodic (1/f) background noise
#'
#' Generates zero-mean noise whose power spectrum falls as `1 / f^exponent` -- the
#' aperiodic ("1/f") background that underlies real electrophysiological recordings
#' and on which oscillatory rhythms sit. It is made by shaping white noise in the
#' frequency domain. Added to a mechanistic oscillator (e.g. [jansenRit()]), it
#' turns an over-concentrated single-rhythm spectrum into a realistic broadband one.
#'
#' @param n Number of samples.
#' @param exponent Spectral slope: power scales as `1 / f^exponent` (1 = pink,
#'   2 = brown; a resting scalp EEG is typically ~1-2).
#' @param sd Target standard deviation of the output.
#' @param seed Optional RNG seed.
#' @return a length-`n` numeric vector.
#' @references He BJ (2014) Scale-free brain activity. Trends Cogn Sci 18:480-487;
#'   Donoghue T et al. (2020) Nat Neurosci 23:1655-1665 (aperiodic/periodic).
#' @seealso [jansenRit()]
#' @export
#' @examples
#' x <- aperiodicNoise(1024, exponent = 1.5, seed = 1)
#' round(sd(x), 3)
aperiodicNoise <- function(n, exponent = 1.5, sd = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  W <- stats::fft(stats::rnorm(n))
  k <- 0:(n - 1); fmag <- pmin(k, n - k)                 # folded frequency magnitude
  amp <- ifelse(fmag == 0, 0, fmag^(-exponent / 2))      # shape to 1/f^exponent power; zero DC
  x <- Re(stats::fft(W * amp, inverse = TRUE)) / n
  x <- x - mean(x); s <- stats::sd(x)
  if (s > 0) x <- x / s * sd
  x
}

#' Jansen-Rit neural-mass EEG generator
#'
#' Simulates a single cortical column with the Jansen-Rit neural-mass model: three
#' interacting neural populations (pyramidal, excitatory and inhibitory
#' interneurons) whose average post-synaptic potentials evolve as second-order
#' ("alpha-function") systems coupled through a sigmoidal rate-to-potential
#' nonlinearity, driven by a background firing-rate input `p(t)`. With the default
#' parameters the column produces a self-organised ~10 Hz **alpha** rhythm; the
#' scalp EEG is taken as the pyramidal post-synaptic potential `y1 - y2`.
#'
#' The excitation/inhibition balance and the drive select the regime: reducing the
#' mean input `p_mean` below the oscillation threshold collapses the alpha limit
#' cycle to a low-amplitude resting fixed point, while changing the gains `A`/`B`
#' or the connectivity `C` shifts the rhythm.
#'
#' @param duration Recording length (s).
#' @param sfeeg Output sampling rate (Hz).
#' @param A,B Excitatory and inhibitory average synaptic gains (mV).
#' @param a,b Excitatory and inhibitory rate constants (inverse time constants, 1/s).
#' @param C Global connectivity constant; the four intra-column connectivities are
#'   `C1 = C`, `C2 = 0.8 C`, `C3 = C4 = 0.25 C`.
#' @param p_mean,p_sd Mean and SD of the background firing-rate input `p(t)`
#'   (pulses/s), drawn independently each integration step. The default
#'   `p_mean = 220` sits in the middle of Jansen & Rit's classic 120-320 input
#'   range, where the column settles into its alpha limit cycle.
#' @param oversample Internal integration steps per output sample (accuracy
#'   margin for the fast `a`-dynamics).
#' @param aperiodic Amplitude of an added 1/f background, as a multiple of the
#'   oscillation's SD (0 = the pure neural-mass rhythm; a real scalp channel needs
#'   `aperiodic` around 1, which spreads power off the alpha peak so the relative
#'   band powers match a real EEG). See [aperiodicNoise()].
#' @param aperiodic_exp Spectral slope of the 1/f background (see [aperiodicNoise()]).
#' @param seed Optional RNG seed.
#' @return a `jansenRit` object: `time`, `eeg` (the `y1 - y2` scalp potential, mV),
#'   `dominant_hz` (the dominant rhythm above 3 Hz, robust to the 1/f floor),
#'   `band` (its conventional EEG band), `sfeeg` and `params`.
#' @references Jansen BH, Rit VG (1995) Electroencephalogram and visual evoked
#'   potential generation in a mathematical model of coupled cortical columns.
#'   Biol Cybern 73:357-366. Grimbert F, Faugeras O (2006) Neural Comput 18:3052-3068.
#' @seealso [ecgsyn()], [kuramoto()], [cpgMatsuoka()]
#' @export
#' @examples
#' eeg <- jansenRit(duration = 6, seed = 1)
#' eeg$dominant_hz          # ~10 Hz alpha rhythm
#' eeg$band                 # "alpha"
jansenRit <- function(duration = 8, sfeeg = 256, A = 3.25, B = 22,
                      a = 100, b = 50, C = 135, p_mean = 220, p_sd = 30,
                      oversample = 4, aperiodic = 0, aperiodic_exp = 1.4, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  C1 <- C; C2 <- 0.8 * C; C3 <- 0.25 * C; C4 <- 0.25 * C
  e0 <- 2.5; v0 <- 6; r <- 0.56
  sig <- function(v) 2 * e0 / (1 + exp(r * (v0 - v)))     # rate <- potential
  dt <- 1 / (sfeeg * oversample)
  n_out <- round(duration * sfeeg)
  n_int <- n_out * oversample
  p <- p_mean + p_sd * stats::rnorm(n_int)                # background drive per step
  # state y = (y0, y1, y2, y0', y1', y2'); scalp EEG = y1 - y2
  deriv <- function(y, pin) {
    c(y[4], y[5], y[6],
      A * a * sig(y[2] - y[3])            - 2 * a * y[4] - a^2 * y[1],
      A * a * (pin + C2 * sig(C1 * y[1])) - 2 * a * y[5] - a^2 * y[2],
      B * b * (C4 * sig(C3 * y[1]))       - 2 * b * y[6] - b^2 * y[3])
  }
  y <- c(0.1, 0.1, 0.1, 0, 0, 0)                          # small non-zero init
  eeg <- numeric(n_out); k <- 0L
  for (i in seq_len(n_int)) {
    pin <- p[i]
    k1 <- deriv(y, pin);          k2 <- deriv(y + dt/2 * k1, pin)
    k3 <- deriv(y + dt/2 * k2, pin); k4 <- deriv(y + dt * k3, pin)
    y <- y + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
    if (i %% oversample == 0L) { k <- k + 1L; eeg[k] <- y[2] - y[3] }
  }
  if (aperiodic > 0)                                      # add a realistic 1/f background
    eeg <- eeg + aperiodicNoise(n_out, aperiodic_exp, sd = aperiodic * stats::sd(eeg))
  t <- seq_len(n_out) / sfeeg
  psd <- .eeg_psd(eeg, sfeeg)
  sel <- psd$freq > 3                                     # the dominant RHYTHM, above the 1/f floor
  dom <- psd$freq[sel][which.max(psd$spec[sel])]
  structure(list(time = t, eeg = eeg, dominant_hz = dom, band = .eeg_band(dom),
                 sfeeg = sfeeg,
                 params = list(A = A, B = B, a = a, b = b, C = C,
                               p_mean = p_mean, p_sd = p_sd,
                               aperiodic = aperiodic, aperiodic_exp = aperiodic_exp)),
            class = "jansenRit")
}

#' @export
print.jansenRit <- function(x, ...) {
  cat(sprintf("Jansen-Rit column -- %.1f s @ %d Hz, dominant %.1f Hz (%s), EEG range [%.2f, %.2f] mV\n",
              max(x$time), x$sfeeg, x$dominant_hz, x$band,
              min(x$eeg), max(x$eeg)))
  invisible(x)
}

#' Relative band power of a simulated EEG
#'
#' Convenience helper: the fraction of spectral power of a `jansenRit` EEG (after
#' dropping the transient) that falls in each conventional band.
#'
#' @param x a `jansenRit` object.
#' @return a named numeric vector of relative power in `delta`, `theta`, `alpha`,
#'   `beta` and `gamma`.
#' @seealso [jansenRit()]
#' @export
#' @examples
#' bandPower(jansenRit(duration = 6, seed = 1))["alpha"]   # alpha-dominant
bandPower <- function(x) {
  stopifnot(inherits(x, "jansenRit"))
  psd <- .eeg_psd(x$eeg, x$sfeeg)
  f <- psd$freq; pw <- psd$spec
  bands <- list(delta = c(0.5, 4), theta = c(4, 8), alpha = c(8, 13),
                beta = c(13, 30), gamma = c(30, Inf))
  tot <- sum(pw)
  vapply(bands, function(bd) sum(pw[f >= bd[1] & f < bd[2]]) / tot, numeric(1))
}
