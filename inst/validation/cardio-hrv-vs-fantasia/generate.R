# cardio-hrv-vs-fantasia: does the closed-loop cardiovascular state-space model,
# fit to a real RR record by WHOLE-RECORD Bayesian inference, reproduce that
# record's heart-rate variability (its variance and its low-/high-frequency
# autonomic balance)? The linear-Gaussian model (a damped Mayer oscillator + a
# respiratory RSA term) is fit to a real Fantasia RR series with
# personalizeCardioWaveform(); a synthetic RR is then generated from the fitted
# model and its HRV compared with the real one.
suppressMessages(library(PhysioTwin))
.args <- commandArgs(FALSE); .f <- sub("^--file=", "", .args[grepl("^--file=", .args)])
HERE <- if (length(.f)) dirname(normalizePath(.f)) else getwd()
OUT <- file.path(HERE, "artifacts")
rr <- read.csv(file.path(OUT, "rr.csv"))$rr_ms
T <- length(rr)

# --- real HRV: power spectrum and LF/HF bands (package internals, used identically
# for real and model so the comparison is like-for-like) --------------------------
cr_psd  <- getFromNamespace(".cr_psd",  "PhysioTwin")
cr_band <- getFromNamespace(".cr_bands", "PhysioTwin")
rb   <- cr_band(cumsum(rr) / 1000, rr)
real_sd    <- sd(rr)
real_lf_hf <- rb$lf_power / rb$hf_power
lf_peak <- rb$lf_peak; hf_peak <- rb$hf_peak          # Mayer (LF) and respiratory (HF) rhythms

# --- fit the closed loop to the whole RR record ----------------------------------
om  <- 2 * pi * lf_peak * mean(rr) / 1000              # Mayer angular freq (rad/beat)
fit <- personalizeCardioWaveform(rr, resp_freq = hf_peak, estimate = c("g", "rho", "A_rsa"),
                                 fixed = list(omega_m = om, rr0 = NULL, q_mayer = 1, r_obs = 9),
                                 n_iter = 3500, seed = 1)
p <- fit$estimates; g <- p[["g"]]; rho <- p[["rho"]]; A_rsa <- p[["A_rsa"]]

# --- synthetic RR from the fitted model; HRV averaged over realisations -----------
A    <- matrix(c(2 * rho * cos(om), 1, -rho^2, 0), 2, 2)
respv <- sin(2 * pi * hf_peak * (0:(T - 1)) * mean(rr) / 1000)
gen_one <- function(seed) {
  set.seed(seed); x <- c(0, 0); s <- numeric(T)
  for (k in seq_len(T)) { x <- A %*% x + c(rnorm(1, 0, 1), 0); s[k] <- mean(rr) + g * x[1] + A_rsa * respv[k] + rnorm(1, 0, 3) }
  s
}
reps <- lapply(1:20, gen_one)
mb   <- vapply(reps, function(s) { b <- cr_band(cumsum(s) / 1000, s); c(b$lf_power, b$hf_power, sd(s)) }, numeric(3))
model_lf_hf <- mean(mb[1, ]) / mean(mb[2, ]); model_sd <- mean(mb[3, ])
syn <- reps[[1]]                                      # a representative synthetic record for the figure

summ <- data.frame(
  metric = c("real_sd_ms", "model_sd_ms", "real_lf_hf", "model_lf_hf",
             "real_lf_peak_hz", "real_hf_peak_hz", "fit_g", "fit_rho", "fit_A_rsa", "n_beats"),
  value = c(round(real_sd, 1), round(model_sd, 1), round(real_lf_hf, 2), round(model_lf_hf, 2),
            round(lf_peak, 3), round(hf_peak, 3), round(g, 2), round(rho, 3), round(A_rsa, 2), T))
write.csv(summ, file.path(OUT, "validation_summary.csv"), row.names = FALSE)

# --- PSDs for the figure ----------------------------------------------------------
pr <- cr_psd(cumsum(rr) / 1000, rr); pm <- cr_psd(cumsum(syn) / 1000, syn)
fm <- pr$freq >= 0 & pr$freq <= 0.5
psd <- data.frame(freq = round(pr$freq[fm], 4),
                  real_psd  = round(pr$spec[fm], 2),
                  model_psd = round(approx(pm$freq, pm$spec, pr$freq[fm], rule = 2)$y, 2))
write.csv(psd, file.path(OUT, "psd.csv"), row.names = FALSE)

ragg::agg_png(file.path(OUT, "figure.png"), width = 1650, height = 700, res = 160)
op <- par(mfrow = c(1, 2), mar = c(4.3, 4.5, 3, 1))
matplot(psd$freq, cbind(psd$real_psd, psd$model_psd), type = "l", lwd = c(3, 2.5), lty = 1,
        col = c("#1B4F72", "#C0392B"), xlim = c(0, 0.5),
        xlab = "frequency (Hz)", ylab = "HRV power spectral density",
        main = sprintf("HRV spectrum: real vs fitted closed loop\n(LF/HF real %.1f, model %.1f)", real_lf_hf, model_lf_hf))
abline(v = c(0.04, 0.15, 0.4), col = "grey70", lty = 3)
legend("topright", c("real (Fantasia f1y03)", "fitted closed-loop model"),
       lwd = c(3, 2.5), col = c("#1B4F72", "#C0392B"), bty = "n", cex = 0.85)
plot(seq_len(200), rr[1:200], type = "l", col = "#1B4F72", lwd = 1.5, ylim = range(rr, syn),
     xlab = "beat", ylab = "RR interval (ms)", main = "RR tachogram: real (blue) vs model (red)")
lines(seq_len(200), syn[1:200], col = "#C0392B", lwd = 1.2)
par(op); grDevices::dev.off()

jsonlite::write_json(list(
  real_data = "Fantasia Database, young healthy subject record f1y03 (ECG beat annotations)",
  authors = "Iyengar N, Peng CK, Morin R, Goldberger AL, Lipsitz LA",
  repository = "PhysioNet", doi = "10.13026/C21596", license = "ODC-BY 1.0",
  url = "https://physionet.org/content/fantasia/",
  model = "closed-loop cardiovascular state-space model, PhysioTwin::cardioStateSpace / personalizeCardioWaveform",
  record_used = sprintf("f1y03, 600-beat resting segment @ 250 Hz; mean RR %.0f ms, SD %.1f ms, LF/HF %.2f (Mayer %.3f Hz, resp %.3f Hz). Closed-loop parameters fit to the whole record by particle-MCMC-family Kalman-likelihood inference: g %.1f, rho %.2f, A_rsa %.1f.",
                        mean(rr), real_sd, real_lf_hf, lf_peak, hf_peak, g, rho, A_rsa)),
  file.path(OUT, "data_sources.json"), auto_unbox = TRUE, pretty = TRUE)

cat(sprintf("DONE real sd=%.1f LF/HF=%.2f (LF %.3f, HF %.3f Hz) | fit g=%.1f rho=%.2f A_rsa=%.1f | model sd=%.1f LF/HF=%.2f\n",
            real_sd, real_lf_hf, lf_peak, hf_peak, g, rho, A_rsa, model_sd, model_lf_hf))
