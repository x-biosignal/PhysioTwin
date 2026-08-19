# ecgsyn-vs-real-ecg: does the mechanistic ECGSYN generator reproduce a REAL ECG?
# The generic model is PERSONALISED (its P/Q/R/S/T event amplitudes fit to the
# real average beat -- the cardiac-twin calibration), then compared to MIT-BIH 100.
suppressMessages(library(PhysioTwin))
.args <- commandArgs(FALSE); .f <- sub("^--file=", "", .args[grepl("^--file=", .args)])
HERE <- if (length(.f)) dirname(normalizePath(.f)) else getwd()   # this case's directory
OUT <- file.path(HERE, "artifacts")
real <- jsonlite::read_json(file.path(OUT, "ecg_real.json"), simplifyVector = TRUE)
fs <- real$fs; N <- 200; pre <- round(0.25 * fs); post <- round(0.45 * fs)
real_beat <- real$avg_beat                             # normalised 0..1, R at idx ~72
b_def <- c(0.25, 0.1, 0.1, 0.1, 0.4)

# one normalised ECGSYN beat template (R-aligned, same window/resample as the real one)
syn_template <- function(a, hr = real$hr) {
  e <- ecgsyn(duration = 4, hr = hr, sfecg = fs, hrv_sd = 0, a = a, seed = 1)
  z <- e$ecg; n <- length(z); th <- 0.6 * max(z); last <- -Inf; rpk <- integer(0)
  for (i in 2:(n - 1)) if (z[i] > th && z[i] >= z[i - 1] && z[i] > z[i + 1] && (i - last) > 0.25 * fs) { rpk <- c(rpk, i); last <- i }
  rpk <- rpk[rpk - pre >= 1 & rpk + post <= n]
  s <- rpk[length(rpk) %/% 2 + 1]; seg <- z[(s - pre):(s + post)]
  v <- approx(seq_along(seg), seg, seq(1, length(seg), length.out = N))$y
  (v - min(v)) / (max(v) - min(v))
}
a_def <- c(1.2, -5, 30, -7.5, 0.75)
default_beat <- syn_template(a_def)
default_cor <- cor(real_beat, default_beat)

# PERSONALISE: fit the 5 event amplitudes so the synthetic beat matches the real one
obj <- function(a) sum((syn_template(a) - real_beat)^2)
fit <- optim(a_def, obj, method = "Nelder-Mead", control = list(maxit = 400, reltol = 1e-7))
a_fit <- fit$par
fitted_beat <- syn_template(a_fit)
fitted_cor <- cor(real_beat, fitted_beat)

# full personalised ECGSYN at the real HR + HRV -> check rate + HRV reproduce
e <- ecgsyn(duration = 120, hr = real$hr, sfecg = fs, hrv_sd = real$sdnn, a = a_fit, seed = 3)
z <- e$ecg; n <- length(z); th <- 0.6 * max(z); last <- -Inf; rpk <- integer(0)
for (i in 2:(n - 1)) if (z[i] > th && z[i] >= z[i - 1] && z[i] > z[i + 1] && (i - last) > 0.25 * fs) { rpk <- c(rpk, i); last <- i }
syn_rr <- diff(rpk) / fs; syn_rr <- syn_rr[syn_rr > 0.3 & syn_rr < 1.5]

summ <- data.frame(
  metric = c("morphology_cor_default", "morphology_cor_personalized", "real_hr", "syn_hr",
             "real_sdnn_ms", "syn_sdnn_ms", "real_n_beats"),
  value = c(round(default_cor, 3), round(fitted_cor, 3), round(real$hr, 1),
            round(60 / mean(syn_rr), 1), round(1000 * real$sdnn, 0),
            round(1000 * sd(syn_rr), 0), real$n_beats))
write.csv(summ, file.path(OUT, "validation_summary.csv"), row.names = FALSE)
write.csv(data.frame(pct_window = seq(0, 1, length.out = N), real = round(real_beat, 4),
                     ecgsyn_default = round(default_beat, 4), ecgsyn_personalized = round(fitted_beat, 4)),
          file.path(OUT, "avg_beat.csv"), row.names = FALSE)

ragg::agg_png(file.path(OUT, "figure.png"), width = 1650, height = 700, res = 160)
op <- par(mfrow = c(1, 2), mar = c(4.3, 4.3, 3, 1)); tt <- seq(-0.25, 0.45, length.out = N)
plot(tt, real_beat, type = "l", lwd = 3, col = "#1B4F72", ylim = c(0, 1),
     xlab = "time from R (s)", ylab = "normalised amplitude",
     main = sprintf("Average beat: real vs ECGSYN\n(default r=%.2f -> personalised r=%.2f)", default_cor, fitted_cor))
lines(tt, default_beat, lwd = 1.5, col = "grey55", lty = 2)
lines(tt, fitted_beat, lwd = 2.5, col = "#C0392B")
legend("topright", c("real (MIT-BIH 100)", "ECGSYN default", "ECGSYN personalised"),
       lwd = c(3, 1.5, 2.5), lty = c(1, 2, 1), col = c("#1B4F72", "grey55", "#C0392B"), bty = "n", cex = 0.85)
strip <- z[round(3 * fs):round(9 * fs)]
plot(seq(0, length(strip) - 1) / fs, strip, type = "l", col = "#C0392B",
     xlab = "time (s)", ylab = "ECG (mV)", main = "Personalised ECGSYN strip (matched HR + HRV)")
par(op); grDevices::dev.off()

jsonlite::write_json(list(
  real_data = "MIT-BIH Arrhythmia Database, record 100 (lead MLII)",
  authors = "Moody GB, Mark RG", repository = "PhysioNet", doi = "10.13026/C2F305",
  license = "ODC-BY 1.0", url = "https://physionet.org/content/mitdb/",
  model = "ECGSYN (McSharry et al. 2003), PhysioTwin::ecgsyn",
  record_used = sprintf("record 100: %d normal beats @ %d Hz; mean HR %.1f bpm, SDNN %.0f ms. ECGSYN event amplitudes personalised to the real average beat, generated at the real HR + HRV.",
                        real$n_beats, fs, real$hr, 1000 * real$sdnn)),
  file.path(OUT, "data_sources.json"), auto_unbox = TRUE, pretty = TRUE)

cat(sprintf("DONE morphology default r=%.3f -> personalised r=%.3f | HR real %.1f vs syn %.1f | SDNN real %.0f vs syn %.0f\n",
            default_cor, fitted_cor, real$hr, 60 / mean(syn_rr), 1000 * real$sdnn, 1000 * sd(syn_rr)))
