# A mechanistic cardiac generator: the ECGSYN dynamical model.
#
# Unlike a Gaussian-template ECG, ECGSYN (McSharry et al. 2003) generates the
# waveform from a 3-state dynamical system whose trajectory circles a limit cycle;
# the ECG (z) is pushed up and down by Gaussian "events" at the P, Q, R, S and T
# angles as the trajectory passes them. Heart-rate variability is injected through
# a realistic RR tachogram (a bimodal LF/HF spectrum), so the generated ECG has
# both physiological MORPHOLOGY and physiological HRV -- exactly what a digital
# twin's cardiac layer needs and what validation against real ECG can test.

# RR tachogram with a bimodal (LF + HF) spectrum -> physiological HRV
.ecg_rrprocess <- function(n, meanrr, sdrr, lfhfratio = 0.5,
                           flo = 0.1, fhi = 0.25, flostd = 0.01, fhistd = 0.01, sfrr = 1) {
  if (n %% 2 == 1) n <- n + 1
  w1 <- 2 * pi * flo; w2 <- 2 * pi * fhi
  c1 <- 2 * pi * flostd; c2 <- 2 * pi * fhistd
  sig2 <- 1; sig1 <- lfhfratio
  df <- sfrr / n; w <- (0:(n - 1)) * 2 * pi * df
  Hw <- sig1 * exp(-(w - w1)^2 / (2 * c1^2)) / sqrt(2 * pi * c1^2) +
        sig2 * exp(-(w - w2)^2 / (2 * c2^2)) / sqrt(2 * pi * c2^2)
  half <- n / 2
  Hw0 <- c(Hw[1:half], Hw[half:1])
  Sw <- (sfrr / 2) * sqrt(Hw0)
  ph <- stats::runif(half - 1) * 2 * pi
  ph <- c(0, ph, 0, -rev(ph))
  x <- (1 / n) * Re(stats::fft(Sw * exp(1i * ph), inverse = TRUE))
  xs <- stats::sd(x); if (xs < 1e-12) xs <- 1
  meanrr + x * (sdrr / xs)
}

#' ECGSYN mechanistic ECG generator (McSharry dynamical model)
#'
#' Generates a realistic ECG from the ECGSYN 3-state limit-cycle model: P/Q/R/S/T
#' waves arise as the trajectory passes Gaussian events, and heart-rate
#' variability is injected via a bimodal RR tachogram. The output has physiological
#' morphology AND HRV.
#'
#' @param duration Length (s).
#' @param hr Mean heart rate (bpm).
#' @param sfecg Sampling frequency (Hz).
#' @param hrv_sd RR-interval standard deviation (s) -- the amount of HRV (0 = fixed
#'   rate).
#' @param lfhfratio LF/HF power ratio of the RR spectrum (autonomic balance).
#' @param anoise Additive measurement-noise SD (mV).
#' @param scale_mv Length-2 physiological output range `c(min, max)` in mV that the
#'   raw model output is affinely mapped onto (McSharry default `c(-0.4, 1.2)`);
#'   `NULL` keeps the raw model units.
#' @param theta,a,b P,Q,R,S,T event angles (rad), amplitudes and widths (McSharry
#'   defaults).
#' @param seed Optional RNG seed.
#' @return an `ecgsyn` object: `time`, `ecg`, the `rr` intervals, `beats` (R-peak
#'   times), `sfecg`.
#' @references McSharry PE, Clifford GD, Tarassenko L, Smith LA (2003) IEEE Trans
#'   Biomed Eng 50:289-294.
#' @seealso [jansenRit()], [windkessel()], [rsaTachogram()], [kuramoto()], [validateTwin()]
#' @export
#' @examples
#' e <- ecgsyn(duration = 6, hr = 60, sfecg = 256, seed = 1)
#' range(e$ecg)          # P/Q/R/S/T morphology
#' length(e$beats)       # ~ 6 beats
ecgsyn <- function(duration = 10, hr = 60, sfecg = 256, hrv_sd = 0.04,
                   lfhfratio = 0.5, anoise = 0, scale_mv = c(-0.4, 1.2),
                   theta = c(-60, -15, 0, 15, 90) * pi / 180,
                   a = c(1.2, -5, 30, -7.5, 0.75),
                   b = c(0.25, 0.1, 0.1, 0.1, 0.4), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- round(duration * sfecg); dt <- 1 / sfecg
  meanrr <- 60 / hr
  nbeats <- ceiling(duration / meanrr) + 4
  rr <- .ecg_rrprocess(nbeats, meanrr, sdrr = hrv_sd, lfhfratio = lfhfratio)
  rr <- pmax(rr, 0.3)                                  # physiological floor
  beat_time <- cumsum(rr)
  t <- (seq_len(n) - 1) * dt
  beat_idx <- findInterval(t, c(0, beat_time)) + 1L
  omega_t <- 2 * pi / rr[pmin(length(rr), pmax(1L, beat_idx))]
  deriv <- function(s, om) {
    xx <- s[1]; yy <- s[2]; zz <- s[3]
    alpha <- 1 - sqrt(xx^2 + yy^2); th <- atan2(yy, xx)
    dth <- ((th - theta + pi) %% (2 * pi)) - pi
    dz <- -sum(a * dth * exp(-dth^2 / (2 * b^2))) - zz
    c(alpha * xx - om * yy, alpha * yy + om * xx, dz)
  }
  Z <- numeric(n); s <- c(1, 0, 0)
  for (i in seq_len(n)) {
    om <- omega_t[i]
    k1 <- deriv(s, om); k2 <- deriv(s + dt/2 * k1, om)
    k3 <- deriv(s + dt/2 * k2, om); k4 <- deriv(s + dt * k3, om)
    s <- s + dt/6 * (k1 + 2*k2 + 2*k3 + k4); Z[i] <- s[3]
  }
  if (!is.null(scale_mv))                              # map raw model units -> physiological mV
    Z <- scale_mv[1] + (scale_mv[2] - scale_mv[1]) * (Z - min(Z)) / (max(Z) - min(Z))
  ecg <- Z + stats::rnorm(n, 0, anoise)
  structure(list(time = t, ecg = ecg, rr = rr, beats = beat_time[beat_time <= duration],
                 sfecg = sfecg, hr = hr), class = "ecgsyn")
}

#' @export
print.ecgsyn <- function(x, ...) {
  cat(sprintf("ECGSYN ECG -- %.1f s @ %d Hz, %d beats (mean HR %.0f bpm, SDNN %.0f ms)\n",
              max(x$time), x$sfecg, length(x$beats), x$hr, 1000 * stats::sd(x$rr)))
  invisible(x)
}
