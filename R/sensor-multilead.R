# Layer 2 -- forward sensor model: a multi-lead ECG derived from one cardiac source.
#
# The heart is a single electrical source; each ECG lead is a different projection
# of it, so leads share the same beat timing but differ in amplitude and shape.
# This maps one clean ECG to a set of named leads, each a gain- and
# morphology-weighted mix of the source and its time-derivative, plus lead-specific
# baseline wander, shared mains interference and independent muscle noise.

# per-lead (projection gain, derivative/morphology weight) for the common leads
.lead_defaults <- function(lead) {
  tbl <- list(I = c(0.9, 0.10), II = c(1.0, 0.00), III = c(0.6, -0.10),
              aVR = c(-0.7, 0.05), aVL = c(0.5, 0.12), aVF = c(0.8, -0.05),
              V1 = c(0.5, 0.30), V2 = c(0.8, 0.25), V3 = c(1.1, 0.12),
              V4 = c(1.2, 0.00), V5 = c(1.0, -0.10), V6 = c(0.8, -0.18))
  if (!is.null(tbl[[lead]])) tbl[[lead]] else c(1.0, 0.0)
}

#' Specify a multi-lead ECG lead set
#'
#' Bundles a set of named ECG leads derived from one cardiac source. Each lead has
#' a projection gain and a morphology (derivative-mixing) weight -- so leads differ
#' in shape, not only amplitude -- plus baseline wander, shared mains interference
#' and muscle noise.
#'
#' @param leads Character vector of lead names (e.g. `"II"`, `"V2"`, `"V5"`).
#' @param gains Optional per-lead projection gains; defaults to standard-ish values.
#' @param baseline_amp,baseline_hz Baseline-wander amplitude (mV) and frequency (Hz).
#' @param powerline_hz,powerline_amp Mains frequency (Hz) and amplitude (mV).
#' @param noise Muscle-noise SD (mV).
#' @param fs Sampling rate (Hz).
#' @return an `ecg_lead_set` object.
#' @seealso [ecgMeasureLeads()], [ecgLead()], [ecgsyn()]
#' @export
#' @examples
#' ecgLeadSet(leads = c("II", "V2", "V5"))
ecgLeadSet <- function(leads = c("II", "V2", "V5"), gains = NULL,
                       baseline_amp = 0.15, baseline_hz = 0.25,
                       powerline_hz = 50, powerline_amp = 0.05, noise = 0.01, fs = 250) {
  dv <- vapply(leads, .lead_defaults, numeric(2))          # row 1 = gain, row 2 = morphology
  structure(list(leads = leads, gains = if (is.null(gains)) dv[1, ] else gains,
                 morph = dv[2, ], baseline_amp = baseline_amp, baseline_hz = baseline_hz,
                 powerline_hz = powerline_hz, powerline_amp = powerline_amp,
                 noise = noise, fs = fs), class = "ecg_lead_set")
}

#' Measure a clean ECG as a multi-lead recording (forward model)
#'
#' Projects one clean ECG onto a set of leads: each lead is a gain- and
#' morphology-weighted mix of the source and its derivative, plus independent
#' baseline wander, shared mains interference and independent muscle noise.
#'
#' @param ecg Clean ECG signal (mV), e.g. from [ecgsyn()].
#' @param sensor An `ecg_lead_set` from [ecgLeadSet()].
#' @param seed Optional RNG seed.
#' @return an `ecg_multilead`: `signals` (`n_time x n_leads`, columns named by lead),
#'   `leads`, `time`, `fs`, sensor.
#' @seealso [ecgLeadSet()], [ecgMeasure()]
#' @export
#' @examples
#' e <- ecgsyn(duration = 10, seed = 1)
#' ml <- ecgMeasureLeads(e$ecg, ecgLeadSet(fs = e$sfecg), seed = 1)
#' dim(ml$signals)
ecgMeasureLeads <- function(ecg, sensor = ecgLeadSet(), seed = NULL) {
  stopifnot(inherits(sensor, "ecg_lead_set"))
  n <- length(ecg); fs <- sensor$fs; t <- (seq_len(n) - 1) / fs
  if (!is.null(seed)) set.seed(seed)
  decg <- c(0, diff(ecg)); s_e <- stats::sd(ecg)
  decg <- decg * s_e / (stats::sd(decg) + 1e-9)             # scale derivative to the ECG's amplitude
  pl <- sensor$powerline_amp * sin(2 * pi * sensor$powerline_hz * t)   # shared mains
  nl <- length(sensor$leads); sig <- matrix(0, n, nl)
  for (l in seq_len(nl)) {
    base <- sensor$baseline_amp * sin(2 * pi * sensor$baseline_hz * t + l)   # per-lead phase
    sig[, l] <- sensor$gains[l] * (ecg + sensor$morph[l] * decg) + base + pl +
                stats::rnorm(n, 0, sensor$noise)
  }
  colnames(sig) <- sensor$leads
  structure(list(signals = sig, leads = sensor$leads, time = t, fs = fs, sensor = sensor),
            class = "ecg_multilead")
}

#' @export
print.ecg_lead_set <- function(x, ...) {
  cat(sprintf("ECG lead set -- leads %s @ %g Hz (baseline %.3g mV, powerline %g Hz)\n",
              paste(x$leads, collapse = ", "), x$fs, x$baseline_amp, x$powerline_hz))
  invisible(x)
}

#' @export
print.ecg_multilead <- function(x, ...) {
  cat(sprintf("Multi-lead ECG -- %d samples x %d leads (%s) @ %g Hz\n",
              nrow(x$signals), length(x$leads), paste(x$leads, collapse = ","), x$fs))
  invisible(x)
}
