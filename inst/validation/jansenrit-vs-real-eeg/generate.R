# jansenrit-vs-real-eeg: does the mechanistic Jansen-Rit neural-mass model
# reproduce the ALPHA RHYTHM of a real eyes-closed resting EEG? The emergent
# alpha peak is compared with a real occipital channel; the model is lightly
# PERSONALISED (its membrane time constants scaled to match the real alpha
# frequency, analogous to the ECGSYN event-amplitude fit).
suppressMessages(library(PhysioTwin))
.args <- commandArgs(FALSE); .f <- sub("^--file=", "", .args[grepl("^--file=", .args)])
HERE <- if (length(.f)) dirname(normalizePath(.f)) else getwd()
OUT <- file.path(HERE, "artifacts")
meta <- jsonlite::read_json(file.path(OUT, "real_meta.json"), simplifyVector = TRUE)
fs <- meta$fs
real <- read.csv(file.path(OUT, "real_eeg.csv"))$eeg_uV

# Welch PSD (Hann, 8 s segments, 50% overlap, density) -- identical for real and model
welch <- function(x, fs, seg = 8) {
  n <- round(seg * fs); if (n > length(x)) n <- length(x)
  step <- n %/% 2; w <- 0.5 - 0.5 * cos(2 * pi * (0:(n - 1)) / (n - 1)); U <- sum(w^2)
  starts <- seq(1, length(x) - n + 1, by = step); P <- numeric(n %/% 2 + 1)
  for (s0 in starts) {
    sg <- (x[s0:(s0 + n - 1)] - mean(x[s0:(s0 + n - 1)])) * w; X <- fft(sg)
    pk <- (Mod(X[1:(n %/% 2 + 1)])^2) / (fs * U); pk[-c(1, length(pk))] <- 2 * pk[-c(1, length(pk))]
    P <- P + pk
  }
  list(freq = (0:(n %/% 2)) * fs / n, psd = P / length(starts))
}
inband <- function(p, lo, hi) p$freq >= lo & p$freq < hi
rel_alpha <- function(p) sum(p$psd[inband(p, 8, 13)]) / sum(p$psd[inband(p, 1, 40)])
dominant  <- function(p) { m <- inband(p, 1, 40); p$freq[m][which.max(p$psd[m])] }
top_band  <- function(p) {                       # which conventional band carries the most power
  b <- c(delta = sum(p$psd[inband(p, 1, 4)]), theta = sum(p$psd[inband(p, 4, 8)]),
         alpha = sum(p$psd[inband(p, 8, 13)]), beta = sum(p$psd[inband(p, 13, 30)]),
         gamma = sum(p$psd[inband(p, 30, 40)])); names(which.max(b))
}

pr <- welch(real, fs)
real_dom <- dominant(pr); real_ra <- rel_alpha(pr); real_top <- top_band(pr)

# generic Jansen-Rit (default parameters) and its emergent alpha
gen <- jansenRit(duration = 60, sfeeg = fs, seed = 1); pg <- welch(gen$eeg, fs)
gen_dom <- dominant(pg)

# PERSONALISE (1): scale the membrane time constants (a, b) so the model's alpha
# peak matches the real one (grid over the alpha regime, dominant averaged over seeds)
mdom <- function(k) mean(vapply(1:3, function(s)
  dominant(welch(jansenRit(duration = 50, sfeeg = fs, a = 100 * k, b = 50 * k, seed = s)$eeg, fs)), numeric(1)))
ks <- seq(0.90, 1.04, by = 0.02); ds <- vapply(ks, mdom, numeric(1))
kbest <- ks[which.min(abs(ds - real_dom))]

# a single cortical column with no broadband background OVER-concentrates alpha:
per0 <- jansenRit(duration = 60, sfeeg = fs, a = 100 * kbest, b = 50 * kbest, seed = 1)
per_ra0 <- rel_alpha(welch(per0$eeg, fs))          # relative alpha with NO 1/f background

# PERSONALISE (2): fit the 1/f (aperiodic) background amplitude so the model's
# RELATIVE ALPHA matches the real one -- a real scalp channel's broadband aperiodic
# background spreads power off the alpha peak (jansenRit(aperiodic=) / aperiodicNoise())
mra <- function(ap) mean(vapply(1:3, function(s)
  rel_alpha(welch(jansenRit(duration = 50, sfeeg = fs, a = 100 * kbest, b = 50 * kbest,
                            aperiodic = ap, seed = s)$eeg, fs)), numeric(1)))
aps <- seq(0, 3, by = 0.2); ras <- vapply(aps, mra, numeric(1))
apbest <- aps[which.min(abs(ras - real_ra))]

per <- jansenRit(duration = 60, sfeeg = fs, a = 100 * kbest, b = 50 * kbest,
                 aperiodic = apbest, seed = 1)
pp <- welch(per$eeg, fs); per_dom <- dominant(pp); per_ra <- rel_alpha(pp); per_top <- top_band(pp)

summ <- data.frame(
  metric = c("real_dominant_hz", "model_dominant_default", "model_dominant_personalized",
             "real_rel_alpha", "model_rel_alpha", "model_rel_alpha_noaperiodic",
             "real_alpha_is_top", "model_alpha_is_top", "ab_scale", "aperiodic_scale", "fs_hz"),
  value = c(round(real_dom, 2), round(gen_dom, 2), round(per_dom, 2),
            round(real_ra, 3), round(per_ra, 3), round(per_ra0, 3),
            as.integer(real_top == "alpha"), as.integer(per_top == "alpha"),
            round(kbest, 3), round(apbest, 2), round(fs, 0)))
write.csv(summ, file.path(OUT, "validation_summary.csv"), row.names = FALSE)

# PSD overlay (1-40 Hz), normalised to peak for shape comparison
fm <- pr$freq >= 1 & pr$freq <= 40; gm <- pp$freq >= 1 & pp$freq <= 40
psd <- data.frame(freq = round(pr$freq[fm], 3),
                  real_psd = round(pr$psd[fm] / max(pr$psd[fm]), 4),
                  model_psd = round(approx(pp$freq[gm], pp$psd[gm] / max(pp$psd[gm]), pr$freq[fm], rule = 2)$y, 4))
write.csv(psd, file.path(OUT, "psd.csv"), row.names = FALSE)

ragg::agg_png(file.path(OUT, "figure.png"), width = 1650, height = 700, res = 160)
op <- par(mfrow = c(1, 2), mar = c(4.3, 4.3, 3, 1))
plot(psd$freq, psd$real_psd, type = "l", lwd = 3, col = "#1B4F72", xlim = c(1, 30),
     xlab = "frequency (Hz)", ylab = "normalised PSD",
     main = sprintf("Alpha rhythm: real EEG vs Jansen-Rit + 1/f\n(peak %.1f/%.1f Hz, rel-alpha %.2f/%.2f)", real_dom, per_dom, real_ra, per_ra))
lines(psd$freq, psd$model_psd, lwd = 2.5, col = "#C0392B")
abline(v = c(8, 13), col = "grey70", lty = 3)
legend("topright", c("real (eegmmidb S001 O1, eyes closed)", "Jansen-Rit (personalised + 1/f)"),
       lwd = c(3, 2.5), col = c("#1B4F72", "#C0392B"), bty = "n", cex = 0.8)
strip <- per$eeg[round(3 * fs):round(8 * fs)]
plot(seq(0, length(strip) - 1) / fs, strip, type = "l", col = "#C0392B",
     xlab = "time (s)", ylab = "EEG (mV)", main = "Personalised Jansen-Rit strip (alpha rhythm)")
par(op); grDevices::dev.off()

jsonlite::write_json(list(
  real_data = "EEG Motor Movement/Imagery Database (eegmmidb), subject S001, run R02 (eyes closed), channel O1",
  authors = "Schalk G, McFarland DJ, Hinterberger T, Birbaumer N, Wolpaw JR",
  repository = "PhysioNet", doi = "10.13026/C28G6P", license = "ODC-BY 1.0",
  url = "https://physionet.org/content/eegmmidb/",
  model = "Jansen-Rit neural-mass model (Jansen & Rit 1995), PhysioTwin::jansenRit",
  record_used = sprintf("S001R02 O1 @ %d Hz, %d samples; real alpha peak %.1f Hz (rel-alpha %.2f). Jansen-Rit membrane time constants scaled x%.3f to match the real alpha frequency, plus a fitted 1/f background of amplitude %.2f (jansenRit aperiodic=) to match the real relative alpha; model rel-alpha %.2f (was %.2f without the 1/f background).",
                        round(fs), meta$n_samples, real_dom, real_ra, kbest, apbest, per_ra, per_ra0)),
  file.path(OUT, "data_sources.json"), auto_unbox = TRUE, pretty = TRUE)

cat(sprintf("DONE real alpha %.2f Hz (rel %.2f) | model %.2f Hz, rel-alpha %.2f -> %.2f via 1/f (aperiodic=%.2f), real %.2f\n",
            real_dom, real_ra, per_dom, per_ra0, per_ra, apbest, real_ra))
