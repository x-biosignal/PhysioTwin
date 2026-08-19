# Layer 4a -- in-silico intervention.
#
# Once a twin is personalised, an intervention is a change to its physical
# parameters -- strengthening raises `strength`, spasticity management lowers
# `stiffness`/`damping`, an orthosis raises `stiffness`. Re-simulating the twin
# with the changed parameters PREDICTS how the movement (and, through the sensor
# layer, the measurable signals) would change -- an experiment run in silico.

# kinematic outcomes of a trajectory
.twin_outcomes <- function(traj) {
  c(rom = diff(range(traj$theta)),
    peak_velocity = max(abs(traj$omega)),
    mean_speed = mean(abs(traj$omega)),
    peak_accel = max(abs(traj$alpha)))
}

#' Predict the effect of an in-silico intervention
#'
#' Applies parameter changes to a (personalised) twin -- multiplicatively via
#' `scale` and/or by absolute values via `set` -- re-simulates, and reports how
#' the movement outcomes change. This is the digital twin used as an experiment:
#' "if this subject's muscle were 30% stronger, how would the movement change?".
#'
#' @param twin A `movement_twin` (typically a personalised one).
#' @param scale Named list of MULTIPLICATIVE changes, e.g. `list(strength = 1.3,
#'   stiffness = 0.7)`.
#' @param set Named list of ABSOLUTE parameter values to set.
#' @param duration,dt,theta0,omega0 Simulation settings (as [simulateTwin()]).
#' @return an `insilico_intervention`: `before`, `after` (trajectories),
#'   `outcomes` (before/after/relative-change table) and the modified `twin`.
#' @seealso [personalizeTwin()], [simulateTwin()]
#' @export
#' @examples
#' tw <- limbTwin(strength = 5, damping = 0.6, stiffness = 3)
#' iv <- insilicoIntervention(tw, scale = list(strength = 1.3, damping = 0.7))
#' iv$outcomes
insilicoIntervention <- function(twin, scale = list(), set = list(),
                                 duration = 6, dt = 0.005, theta0 = 0.1, omega0 = 0) {
  stopifnot(inherits(twin, "movement_twin"))
  post <- twin
  for (nm in names(scale)) post[[nm]] <- post[[nm]] * scale[[nm]]
  for (nm in names(set))   post[[nm]] <- set[[nm]]
  before <- simulateTwin(twin, duration, dt, theta0, omega0)
  after  <- simulateTwin(post, duration, dt, theta0, omega0)
  ob <- .twin_outcomes(before); oa <- .twin_outcomes(after)
  outcomes <- data.frame(outcome = names(ob), before = round(ob, 4),
                         after = round(oa, 4),
                         pct_change = round(100 * (oa - ob) / abs(ob), 1),
                         row.names = NULL)
  structure(list(before = before, after = after, outcomes = outcomes, twin = post,
                 changes = c(scale = list(scale), set = list(set))),
            class = "insilico_intervention")
}

#' @export
print.insilico_intervention <- function(x, ...) {
  cat("In-silico intervention -- predicted movement change\n")
  for (i in seq_len(nrow(x$outcomes)))
    cat(sprintf("  %-14s %8.3f -> %8.3f  (%+.1f%%)\n", x$outcomes$outcome[i],
                x$outcomes$before[i], x$outcomes$after[i], x$outcomes$pct_change[i]))
  invisible(x)
}
