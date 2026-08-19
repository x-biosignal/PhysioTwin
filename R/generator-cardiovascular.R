# Cardiovascular generators: the Windkessel arterial load and a closed-loop
# baroreflex.
#
# Arterial blood pressure is not a template either; it is the response of a
# compliant, resistive arterial tree to the pulsatile flow ejected by the heart.
#   * windkessel()  -- the classic 2-/3-element Windkessel: the arteries as a
#                      compliance C draining through a peripheral resistance R
#                      (3-element adds a proximal characteristic impedance Zc).
#                      A half-sine ejection each beat charges the compliance;
#                      during diastole the pressure decays exponentially with
#                      time constant RC. Produces a realistic pressure wave with
#                      the right systolic/diastolic/mean/pulse pressures.
#   * baroreflex()  -- the short-term pressure-regulating reflex closed around the
#                      Windkessel: each beat the mean pressure is compared to a
#                      set-point and the next RR interval (heart period) is
#                      adjusted (higher pressure -> longer RR, slower heart), so a
#                      pressure perturbation is actively buffered.
# Dependency-free base R, RK4 integration.

# Half-sine aortic ejection flow (mL/s) at within-beat time `ph`, for a beat of
# period `T` with ejection duration `Ts` delivering stroke volume `sv`.
.wk_flow <- function(ph, Ts, sv) {
  imax <- sv * pi / (2 * Ts)                 # so integral over [0,Ts] == sv
  ifelse(ph < Ts, imax * sin(pi * ph / Ts), 0)
}

# Integrate one beat of the Windkessel compliance node forward from pressure
# `pc0`, returning the per-sample node pressure and flow. dPc/dt = (I - Pc/R)/C.
.wk_beat <- function(pc0, T, Ts, sv, R, C, dt) {
  nb <- max(1L, round(T / dt))
  pc <- pc0; Pc <- numeric(nb); Fl <- numeric(nb)
  for (j in seq_len(nb)) {
    tt <- (j - 1) * dt
    f  <- function(p, ph) (.wk_flow(ph, Ts, sv) - p / R) / C
    k1 <- f(pc,             tt)
    k2 <- f(pc + dt/2 * k1, tt + dt/2)
    k3 <- f(pc + dt/2 * k2, tt + dt/2)
    k4 <- f(pc + dt   * k3, tt + dt)
    pc <- pc + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
    Pc[j] <- pc; Fl[j] <- .wk_flow(tt, Ts, sv)
  }
  list(pc = Pc, flow = Fl, end = pc)
}

#' Windkessel arterial blood-pressure generator
#'
#' Simulates arterial pressure with the classic Windkessel model: a half-sine
#' ventricular ejection each beat charges the arterial compliance `C`, which
#' drains through the peripheral resistance `R`; the 3-element variant adds a
#' proximal characteristic impedance `Zc`. During diastole (no inflow) the
#' pressure decays exponentially with time constant `R * C`. With physiological
#' parameters it reproduces the systolic/diastolic/mean/pulse pressures and the
#' diastolic run-off.
#'
#' @param duration Recording length (s).
#' @param hr Heart rate (bpm).
#' @param stroke_volume Stroke volume ejected each beat (mL).
#' @param R Peripheral (systemic vascular) resistance (mmHg*s/mL).
#' @param C Arterial compliance (mL/mmHg).
#' @param Zc Proximal characteristic impedance (mmHg*s/mL); used when `elements = 3`.
#' @param systolic_frac Fraction of the cardiac cycle occupied by ejection.
#' @param elements 2 or 3 (2-element = pure RC; 3-element adds `Zc`).
#' @param sfp Output sampling rate (Hz).
#' @param p0 Initial compliance-node pressure (mmHg).
#' @return a `windkessel` object: `time`, `pressure` (arterial pressure, mmHg),
#'   `flow` (aortic inflow, mL/s), and the derived `systolic`, `diastolic`, `map`
#'   (mean arterial pressure) and `pp` (pulse pressure), measured over the last
#'   settled beats.
#' @references Westerhof N, Lankhaar J-W, Westerhof BE (2009) The arterial
#'   Windkessel. Med Biol Eng Comput 47:131-141.
#' @seealso [baroreflex()], [ecgsyn()], [respiration()]
#' @export
#' @examples
#' bp <- windkessel(duration = 10, hr = 70)
#' round(c(bp$systolic, bp$diastolic))     # ~ 120 / 80 mmHg
windkessel <- function(duration = 10, hr = 70, stroke_volume = 80,
                       R = 1.0, C = 1.9, Zc = 0.05, systolic_frac = 0.3,
                       elements = 3, sfp = 250, p0 = 80) {
  stopifnot(elements %in% c(2, 3))
  dt <- 1 / sfp; T <- 60 / hr; Ts <- systolic_frac * T
  n_beats <- ceiling(duration / T) + 1               # +1: per-beat rounding never undershoots n
  Pc <- numeric(0); Fl <- numeric(0); pc <- p0
  for (b in seq_len(n_beats)) {
    beat <- .wk_beat(pc, T, Ts, stroke_volume, R, C, dt)
    Pc <- c(Pc, beat$pc); Fl <- c(Fl, beat$flow); pc <- beat$end
  }
  n <- round(duration * sfp); Pc <- Pc[seq_len(n)]; Fl <- Fl[seq_len(n)]
  pressure <- if (elements == 3) Pc + Zc * Fl else Pc
  t <- seq_len(n) * dt
  settled <- t > (duration - 3 * T)                 # last ~3 beats, past transient
  structure(list(time = t, pressure = pressure, flow = Fl,
                 systolic = max(pressure[settled]), diastolic = min(pressure[settled]),
                 map = mean(pressure[settled]),
                 pp = max(pressure[settled]) - min(pressure[settled]),
                 hr = hr, R = R, C = C, tau = R * C),
            class = "windkessel")
}

#' @export
print.windkessel <- function(x, ...) {
  cat(sprintf("Windkessel -- %.0f/%.0f mmHg (MAP %.0f, PP %.0f), HR %.0f, tau=RC %.2f s\n",
              x$systolic, x$diastolic, x$map, x$pp, x$hr, x$tau))
  invisible(x)
}

#' Baroreflex closed-loop blood-pressure regulation
#'
#' Closes the short-term baroreflex around the [windkessel()] arterial model:
#' each beat the mean arterial pressure is compared to a set-point and the next
#' heart period (RR interval) is set to `RR0 + gain * (MAP - setpoint)` -- a
#' higher pressure lengthens RR (slows the heart), lowering cardiac output and
#' buffering the pressure (negative feedback; `gain` is the baroreflex
#' sensitivity, s/mmHg). A step change in peripheral resistance at `perturb_time`
#' (e.g. vasodilation) challenges the loop; with the reflex engaged the mean
#' pressure is actively regulated back toward the set-point.
#'
#' @param duration Recording length (s).
#' @param setpoint Regulated mean-arterial-pressure set-point (mmHg).
#' @param gain Baroreflex sensitivity, RR change per unit pressure error (s/mmHg).
#' @param hr0 Baseline heart rate (bpm).
#' @param stroke_volume,R,C,Zc,systolic_frac Windkessel parameters (see [windkessel()]).
#' @param perturb_time Time of the resistance step (s); `Inf` for no perturbation.
#' @param perturb_R Multiplier applied to `R` from `perturb_time` on (e.g. 0.7 =
#'   30% vasodilation).
#' @param reflex If `FALSE`, the RR interval is held at baseline (open loop) -- the
#'   control for demonstrating regulation.
#' @param sfp Output sampling rate (Hz).
#' @return a `baroreflex` object: continuous `time`/`pressure`; per-beat
#'   `beat_time`, `rr`, `hr`, `beat_map`; the `setpoint`; and `regulated` (the
#'   post-perturbation steady-state absolute MAP error).
#' @references Ursino M (1998) Am J Physiol 275:H1733-H1747. deBoer RW et al.
#'   (1987) Am J Physiol 253:H680-H689.
#' @seealso [windkessel()], [rsaTachogram()]
#' @export
#' @examples
#' br <- baroreflex(duration = 60, perturb_time = 30, perturb_R = 0.7)
#' br$regulated       # small residual MAP error -- pressure is buffered
baroreflex <- function(duration = 60, setpoint = 92, gain = 0.02, hr0 = 70,
                       stroke_volume = 80, R = 1.0, C = 1.9, Zc = 0.05,
                       systolic_frac = 0.3, perturb_time = 30, perturb_R = 0.7,
                       reflex = TRUE, sfp = 250) {
  dt <- 1 / sfp; RR0 <- 60 / hr0; RR <- RR0
  pc <- setpoint; tg <- 0
  P <- numeric(0); bt <- numeric(0); rrs <- numeric(0); maps <- numeric(0)
  while (tg < duration) {
    Rcur <- if (tg >= perturb_time) R * perturb_R else R
    T <- RR; Ts <- systolic_frac * T
    beat <- .wk_beat(pc, T, Ts, stroke_volume, Rcur, C, dt)
    press <- beat$pc + Zc * beat$flow; pc <- beat$end
    mapb <- mean(press)
    P <- c(P, press); bt <- c(bt, tg); rrs <- c(rrs, RR); maps <- c(maps, mapb)
    tg <- tg + T
    if (reflex) RR <- max(0.3, min(2.0, RR0 + gain * (mapb - setpoint)))  # BR update
  }
  n <- min(round(duration * sfp), length(P)); P <- P[seq_len(n)]
  post <- bt >= perturb_time & bt <= duration          # after the challenge
  mp   <- maps[post]                                    # settled-error over last 5 post beats
  reg  <- if (length(mp)) abs(mean(mp[max(1L, length(mp) - 4L):length(mp)]) - setpoint) else NA_real_
  structure(list(time = seq_len(n) * dt, pressure = P, beat_time = bt,
                 rr = rrs, hr = 60 / rrs, beat_map = maps, setpoint = setpoint,
                 regulated = reg, reflex = reflex),
            class = "baroreflex")
}

#' @export
print.baroreflex <- function(x, ...) {
  cat(sprintf("Baroreflex (%s) -- setpoint %.0f mmHg, HR %.0f-%.0f bpm, residual MAP error %.1f mmHg\n",
              if (x$reflex) "closed loop" else "open loop", x$setpoint,
              min(x$hr), max(x$hr), x$regulated))
  invisible(x)
}
