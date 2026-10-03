# Respiration generator and respiratory sinus arrhythmia (RSA).
#
# Breathing and its coupling to the heartbeat are mechanistic, not templates:
#   * respiration()   -- an asymmetric breathing waveform (inspiration shorter
#                        than expiration, set by the I:E ratio); chest volume
#                        rises on inspiration and falls on expiration.
#   * rsaTachogram()  -- respiratory sinus arrhythmia by INTEGRAL PULSE FREQUENCY
#                        MODULATION (the standard mechanistic model of heartbeat
#                        timing): the autonomic outflow modulates the instantaneous
#                        heart rate -- faster on inspiration (RSA), plus a ~0.1 Hz
#                        Mayer wave -- and a beat is emitted each time the integral
#                        of that rate crosses an integer. The resulting RR interval
#                        series (tachogram) carries the high-frequency (respiratory)
#                        and low-frequency (Mayer) peaks of heart-rate variability.
# Dependency-free base R.

#' Respiratory waveform generator
#'
#' Simulates a breathing (chest-volume) waveform: within each breath the volume
#' rises smoothly on inspiration and falls on expiration, with the inspiration
#' fraction set by the inspiration:expiration (`ie_ratio`) ratio (physiological
#' breathing is expiration-dominant, ~1:2).
#'
#' @param duration Recording length (s).
#' @param rate Breathing rate (breaths/min).
#' @param ie_ratio Inspiration-to-expiration duration ratio (0.5 = 1:2).
#' @param amplitude Peak-to-trough volume amplitude (arbitrary units).
#' @param sfp Sampling rate (Hz).
#' @return a `respiration` object: `time`, `volume` (chest expansion, zero-mean),
#'   `flow` (its time derivative), and `rate`.
#' @seealso [rsaTachogram()], [windkessel()]
#' @export
#' @examples
#' rp <- respiration(duration = 30, rate = 15)
#' rp$rate
respiration <- function(duration = 60, rate = 15, ie_ratio = 0.5,
                        amplitude = 1, sfp = 100) {
  dt <- 1 / sfp; n <- round(duration * sfp); t <- (seq_len(n) - 1) * dt
  Tb <- 60 / rate                                   # breath period
  fi <- ie_ratio / (1 + ie_ratio)                   # inspiration fraction of the breath
  Ti <- fi * Tb; Te <- Tb - Ti
  ph <- t %% Tb
  vol <- ifelse(ph < Ti,
                0.5 * (1 - cos(pi * ph / Ti)),                 # inspiration 0 -> 1
                0.5 * (1 + cos(pi * (ph - Ti) / Te)))          # expiration 1 -> 0
  vol <- amplitude * (vol - mean(vol))                         # zero-mean
  flow <- c(0, diff(vol)) / dt
  structure(list(time = t, volume = vol, flow = flow, rate = rate,
                 ie_ratio = ie_ratio, sfp = sfp),
            class = "respiration")
}

#' @export
print.respiration <- function(x, ...) {
  cat(sprintf("Respiration -- %.0f breaths/min, I:E 1:%.1f, %.1f s @ %d Hz\n",
              x$rate, 1 / x$ie_ratio, max(x$time), x$sfp))
  invisible(x)
}

# Smoothed HRV power spectrum of an RR tachogram: interpolate RR onto a uniform
# grid, detrend, and return a smoothed one-sided PSD (Hz vs power).
.hrv_psd <- function(beat_time, rr, fs_i = 4) {
  tt <- beat_time[-1]                                # RR[i] is dated at the i-th beat
  grid <- seq(min(tt), max(tt), by = 1 / fs_i)
  ri <- stats::approx(tt, rr, xout = grid, rule = 2)$y
  ri <- ri - mean(ri)
  sp <- stats::spec.pgram(ri, spans = c(5, 5), taper = 0.1, detrend = TRUE, plot = FALSE)
  list(freq = sp$freq * fs_i, spec = sp$spec)
}

#' Respiratory sinus arrhythmia tachogram (heart-rate variability generator)
#'
#' Generates a heart-rate-variability tachogram by integral pulse frequency
#' modulation: the instantaneous heart rate is modulated by respiration (faster
#' on inspiration -- respiratory sinus arrhythmia) and by a ~0.1 Hz Mayer wave,
#' and a beat is emitted whenever the integral of the instantaneous rate crosses
#' an integer. The resulting RR-interval series reproduces the high-frequency
#' (respiratory) and low-frequency (Mayer) peaks of the HRV spectrum.
#'
#' @param duration Recording length (s).
#' @param hr0 Mean heart rate (bpm).
#' @param resp_rate Breathing rate (breaths/min); sets the HF peak.
#' @param rsa Depth of the respiratory heart-rate modulation (fraction of `hr0`).
#' @param mayer Depth of the Mayer-wave modulation (fraction of `hr0`).
#' @param mayer_freq Mayer-wave frequency (Hz).
#' @param ie_ratio Inspiration:expiration ratio of the driving respiration.
#' @param noise SD of a white heart-rate jitter (fraction of `hr0`).
#' @param seed Optional RNG seed.
#' @return an `rsaTachogram` object: `beat_time`, `rr` (RR intervals, s), `hr`
#'   (instantaneous, bpm), the driving `resp` waveform; and the HRV summaries
#'   `hf_peak` (Hz of the HF spectral peak), `lf_power`, `hf_power` and `lf_hf`
#'   (the LF/HF ratio).
#' @references Berger RD et al. (1989) IEEE Trans Biomed Eng 36:1049-1055
#'   (IPFM). Task Force (1996) Circulation 93:1043-1065 (HRV bands).
#' @seealso [respiration()], [ecgsyn()]
#' @export
#' @examples
#' hrv <- rsaTachogram(duration = 300, resp_rate = 15, seed = 1)
#' hrv$hf_peak            # ~ 0.25 Hz, the respiratory rate
rsaTachogram <- function(duration = 300, hr0 = 70, resp_rate = 15, rsa = 0.12,
                         mayer = 0.06, mayer_freq = 0.1, ie_ratio = 0.5,
                         noise = 0.01, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  fs <- 100; dt <- 1 / fs; n <- round(duration * fs); t <- (seq_len(n) - 1) * dt
  rn <- respiration(duration, resp_rate, ie_ratio, sfp = fs)$volume
  rn <- rn / max(abs(rn))                                        # unit-amplitude modulator
  f0 <- hr0 / 60                                                 # mean rate (Hz)
  drive <- f0 * (1 + rsa * rn + mayer * sin(2 * pi * mayer_freq * t))
  if (noise > 0) drive <- drive + f0 * noise * stats::rnorm(n)
  drive <- pmax(drive, 1e-3)                                     # rate stays positive
  phase <- cumsum(drive) * dt                                    # IPFM integral
  beat_idx <- which(diff(floor(phase)) >= 1)                     # integer crossings = beats
  beat_time <- t[beat_idx]
  rr <- diff(beat_time)
  hf <- .hrv_psd(beat_time, rr)
  inb <- function(lo, hi) hf$freq >= lo & hf$freq < hi
  lf_p <- sum(hf$spec[inb(0.04, 0.15)]); hf_p <- sum(hf$spec[inb(0.15, 0.40)])
  hf_peak <- hf$freq[inb(0.15, 0.40)][which.max(hf$spec[inb(0.15, 0.40)])]
  structure(list(beat_time = beat_time, rr = rr, hr = 60 / rr,
                 resp = rn, resp_rate = resp_rate, time = t,
                 hf_peak = hf_peak, lf_power = lf_p, hf_power = hf_p,
                 lf_hf = lf_p / hf_p),
            class = "rsaTachogram")
}

#' @export
print.rsaTachogram <- function(x, ...) {
  cat(sprintf("RSA tachogram -- %d beats, mean HR %.0f bpm, HF peak %.2f Hz (resp %.2f Hz), LF/HF %.2f\n",
              length(x$rr), mean(x$hr), x$hf_peak, x$resp_rate / 60, x$lf_hf))
  invisible(x)
}
