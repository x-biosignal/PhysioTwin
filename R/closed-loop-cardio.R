# A closed-loop cardiovascular-respiratory twin.
#
# The earlier generators are open loop: windkessel() takes a fixed heart rate,
# baroreflex() controls only heart rate, and rsaTachogram() IMPOSES the
# respiratory and Mayer rhythms on the heartbeat. Real cardiovascular variability
# is EMERGENT: the baroreflex regulates both heart period and peripheral
# resistance with a delay, and that delayed negative feedback loop, excited by
# noise, spontaneously oscillates near 0.1 Hz -- the Mayer wave -- with no
# oscillator imposed. Respiration, coupled in mechanically (stroke volume) and
# centrally (vagal tone), adds the respiratory-frequency oscillations of blood
# pressure (Traube-Hering waves) and heart rate (respiratory sinus arrhythmia).
#
# A DeBoer-style beat-to-beat model: a Windkessel pressure map closed by a
# two-arm (heart-period + resistance), delayed baroreflex, driven by respiration.
# Dependency-free base R.

# smoothed one-sided PSD of a beat-series `x` dated at `bt`, resampled at fs Hz
.cr_psd <- function(bt, x, fs = 4) {
  grid <- seq(bt[1], bt[length(bt)], by = 1 / fs)
  xi <- stats::approx(bt, x, xout = grid, rule = 2)$y; xi <- xi - mean(xi)
  sp <- stats::spec.pgram(xi, spans = c(7, 7), taper = 0.1, detrend = TRUE, plot = FALSE)
  list(freq = sp$freq * fs, spec = sp$spec)
}

# LF/HF band summary of a beat series (LF = Mayer band, HF = respiratory band)
.cr_bands <- function(bt, x) {
  p <- .cr_psd(bt, x); f <- p$freq
  lf <- f >= 0.04 & f < 0.15; hf <- f >= 0.15 & f < 0.40
  list(lf_power = sum(p$spec[lf]), hf_power = sum(p$spec[hf]),
       lf_peak = f[lf][which.max(p$spec[lf])], hf_peak = f[hf][which.max(p$spec[hf])])
}

#' Closed-loop cardiovascular-respiratory twin
#'
#' A beat-to-beat model in which arterial pressure, heart period and peripheral
#' resistance evolve under a delayed, two-arm baroreflex, coupled to respiration.
#' A Windkessel map sets each beat's systolic and diastolic pressure; the
#' baroreflex adjusts the next heart period (fast vagal + delayed sympathetic) and
#' the peripheral resistance (delayed sympathetic) from the sensed systolic
#' pressure; respiration modulates stroke volume (a mechanical blood-pressure
#' effect) and vagal tone (respiratory sinus arrhythmia). The delayed
#' resistance/heart-period feedback, excited by small noise, produces an EMERGENT
#' ~0.1 Hz Mayer wave in pressure and heart rate; respiration adds the
#' respiratory-frequency oscillations. None of these rhythms is imposed.
#'
#' @param duration Recording length (s).
#' @param hr0 Baseline heart rate (bpm).
#' @param resp_rate Breathing rate (breaths/min).
#' @param S0,D0 Target systolic/diastolic pressure (mmHg).
#' @param R0 Baseline peripheral resistance (relative).
#' @param g_vagal,g_symp Baroreflex heart-period gains (s/mmHg): fast vagal (current
#'   beat) and delayed sympathetic.
#' @param symp_delay,resist_delay Sympathetic heart-period and resistance baroreflex
#'   delays (beats).
#' @param g_resist Baroreflex resistance gain (per mmHg).
#' @param rsa_amp Respiratory vagal modulation of heart period (s).
#' @param resp_bp_amp Respiratory mechanical modulation of stroke volume (mmHg).
#' @param starling Frank-Starling sensitivity of stroke volume to the preceding
#'   heart period (fractional).
#' @param noise SD of the beat-to-beat neural noise (fraction) that excites the loop.
#' @param respiration If `FALSE`, respiration is switched off (to show the Mayer
#'   wave is intrinsic to the baroreflex loop, not driven by breathing).
#' @param seed Optional RNG seed.
#' @return a `cardio_respiratory` object: per-beat `time`, `rr`, `hr`, `systolic`,
#'   `diastolic`, `resistance`; the HRV and systolic-pressure LF/HF band summaries
#'   `hrv` and `bpv` (each with `lf_power`/`hf_power`/`lf_peak`/`hf_peak`); and the
#'   baroreflex sensitivity `brs` (ms/mmHg) recovered from the beat series.
#' @references de Boer RW, Karemaker JM, Strackee J (1987) Am J Physiol
#'   253:H680-H689; Julien C (2006) Cardiovasc Res 70:12-21 (Mayer waves).
#' @seealso [windkessel()], [baroreflex()], [rsaTachogram()]
#' @export
#' @examples
#' cr <- cardioRespiratory(duration = 240, resp_rate = 12, seed = 1)
#' cr$bpv$lf_peak     # ~0.1 Hz Mayer wave in blood pressure
#' cr$hrv$hf_peak     # respiratory-frequency RSA in heart rate
cardioRespiratory <- function(duration = 300, hr0 = 70, resp_rate = 12,
                              S0 = 120, D0 = 80, R0 = 1,
                              g_vagal = 0.010, g_symp = 0.007, symp_delay = 2L,
                              g_resist = 0.010, resist_delay = 4L,
                              rsa_amp = 0.035, resp_bp_amp = 2.5,
                              starling = 0.10, noise = 0.02,
                              respiration = TRUE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  I0 <- 60 / hr0; PP0 <- S0 - D0; f_r <- resp_rate / 60
  tau <- -I0 / (R0 * log(D0 / S0))              # Windkessel decay so D0 = S0 exp(-I0/(R0 tau))
  nmax <- ceiling(duration / I0 * 1.6) + 20L
  I <- S <- D <- R <- tt <- numeric(nmax)
  I[1:5] <- I0; S[1:5] <- S0; D[1:5] <- D0; R[1:5] <- R0
  for (k in 2:5) tt[k] <- tt[k - 1] + I0
  n <- 5L
  while (tt[n] < duration && n < nmax) {
    n <- n + 1L
    tt[n] <- tt[n - 1] + I[n - 1]
    rsig <- if (respiration) sin(2 * pi * f_r * tt[n]) else 0
    sv <- PP0 * (1 + starling * (I[n - 1] - I0) / I0) + resp_bp_amp * rsig   # Starling + respiration
    S[n] <- D[n - 1] + sv                                                    # systolic builds on diastole
    si <- max(1L, n - symp_delay); ri <- max(1L, n - resist_delay)
    I[n] <- I0 + g_vagal * (S[n] - S0) + g_symp * (S[si] - S0) +             # baroreflex heart period
            rsa_amp * rsig + noise * I0 * stats::rnorm(1)
    I[n] <- max(0.3, min(2.0, I[n]))
    R[n] <- max(0.4, R0 - g_resist * (S[ri] - S0) + noise * stats::rnorm(1)) # baroreflex resistance (delayed)
    D[n] <- S[n] * exp(-I[n] / (R[n] * tau))                                 # diastolic runoff
  }
  k <- 6:n; bt <- tt[k]
  hrv <- .cr_bands(bt, I[k]); bpv <- .cr_bands(bt, S[k])
  brs <- stats::coef(stats::lm(I[k][-1] ~ S[k][-length(k)]))[2] * 1000       # ms/mmHg (I vs prior S)
  structure(list(time = bt, rr = I[k], hr = 60 / I[k], systolic = S[k],
                 diastolic = D[k], resistance = R[k], hrv = hrv, bpv = bpv,
                 brs = unname(brs), respiration = respiration),
            class = "cardio_respiratory")
}

#' @export
print.cardio_respiratory <- function(x, ...) {
  cat(sprintf("Closed-loop cardio-respiratory -- %d beats, HR %.0f, BP %.0f/%.0f mmHg\n",
              length(x$rr), mean(x$hr), mean(x$systolic), mean(x$diastolic)))
  cat(sprintf("  BP Mayer (LF) peak %.3f Hz; HR RSA (HF) peak %.3f Hz; BRS %.1f ms/mmHg%s\n",
              x$bpv$lf_peak, x$hrv$hf_peak, x$brs,
              if (!x$respiration) " (respiration off)" else ""))
  invisible(x)
}
