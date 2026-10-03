# Layer 4c -- multi-modal coupling: one latent driver, many modalities.
#
# Movement and physiology are not independent: a subject's fitness raises movement
# capacity AND lowers resting heart rate, and effort raises both movement demand
# and cardiorespiratory drive. A multi-modal twin couples a latent fitness/effort
# state to BOTH a movement twin and physiological set-points, so an intervention
# that changes fitness produces a physiologically COHERENT change across modalities
# -- e.g. a rehab programme that raises movement range of motion while lowering
# resting heart rate.

#' Build a multi-modal (movement + cardiorespiratory) twin
#'
#' Couples a latent `fitness` and `effort` to both a [limbTwin()] movement model
#' (higher fitness -> more strength, less damping) and cardiorespiratory set-points
#' (higher fitness -> lower resting heart rate; higher effort -> higher heart rate
#' and breathing rate).
#'
#' @param fitness Latent fitness (0..1).
#' @param effort Latent effort/activity level (0..1).
#' @param base_strength,base_damping Reference movement parameters at `fitness = 0.5`.
#' @return a `multimodal_twin`: the derived `movement` twin, `resting_hr`,
#'   `resp_rate` and the latent `fitness`/`effort`.
#' @seealso [simulateMultimodal()], [insilicoInterventionMultimodal()], [rsaTachogram()]
#' @export
#' @examples
#' mt <- multimodalTwin(fitness = 0.7)
#' mt$resting_hr
multimodalTwin <- function(fitness = 0.5, effort = 0.5, base_strength = 5,
                           base_damping = 0.5) {
  fitness <- min(1, max(0, fitness)); effort <- min(1, max(0, effort))
  movement <- limbTwin(strength = base_strength * (0.6 + 0.8 * fitness),
                       damping  = base_damping  * (1.3 - 0.6 * fitness))
  structure(list(movement = movement,
                 resting_hr = 80 - 30 * fitness + 20 * effort,   # fitter = lower resting HR
                 resp_rate = 12 + 12 * effort, fitness = fitness, effort = effort),
            class = "multimodal_twin")
}

#' Simulate a multi-modal twin (movement + cardiorespiratory summary)
#'
#' Runs the movement twin and the cardiorespiratory generators and returns a small
#' cross-modal summary: movement range of motion and peak velocity, and the mean
#' heart rate.
#'
#' @param mt a `multimodal_twin` from [multimodalTwin()].
#' @param duration_move Movement simulation length (s).
#' @param duration_hr Heart-rate simulation length (s).
#' @param dt,theta0 Movement simulation settings.
#' @param seed Optional RNG seed.
#' @return a named numeric vector: `rom`, `peak_velocity`, `mean_hr`, `resp_rate`.
#' @seealso [multimodalTwin()], [insilicoInterventionMultimodal()]
#' @export
#' @examples
#' simulateMultimodal(multimodalTwin(fitness = 0.6), seed = 1)
simulateMultimodal <- function(mt, duration_move = 6, duration_hr = 120,
                               dt = 0.005, theta0 = 0.1, seed = NULL) {
  stopifnot(inherits(mt, "multimodal_twin"))
  tr <- simulateTwin(mt$movement, duration = duration_move, dt = dt, theta0 = theta0)
  o <- .twin_outcomes(tr)
  hrv <- rsaTachogram(duration = duration_hr, hr0 = mt$resting_hr,
                      resp_rate = mt$resp_rate, seed = seed)
  c(rom = unname(o["rom"]), peak_velocity = unname(o["peak_velocity"]),
    mean_hr = mean(hrv$hr), resp_rate = mt$resp_rate)
}

#' Multi-modal in-silico intervention
#'
#' Applies a fitness-raising intervention (e.g. a rehabilitation / conditioning
#' programme) to a multi-modal twin and predicts the coherent change across
#' modalities: movement range of motion up, resting/mean heart rate down.
#'
#' @param mt a `multimodal_twin` from [multimodalTwin()].
#' @param fitness_gain Increase in latent fitness (0..1, clamped).
#' @param ... Passed to [simulateMultimodal()].
#' @return an `mm_intervention`: the `before`/`after` cross-modal summaries and
#'   their `change`.
#' @seealso [multimodalTwin()], [populationIntervention()]
#' @export
#' @examples
#' insilicoInterventionMultimodal(multimodalTwin(fitness = 0.4), fitness_gain = 0.3, seed = 1)$change
insilicoInterventionMultimodal <- function(mt, fitness_gain = 0.2, ...) {
  stopifnot(inherits(mt, "multimodal_twin"))
  post <- multimodalTwin(fitness = mt$fitness + fitness_gain, effort = mt$effort)
  before <- simulateMultimodal(mt, ...); after <- simulateMultimodal(post, ...)
  structure(list(before = before, after = after, change = after - before,
                 fitness_gain = fitness_gain), class = "mm_intervention")
}

#' @export
print.multimodal_twin <- function(x, ...) {
  cat(sprintf("Multi-modal twin -- fitness %.2f, effort %.2f -> strength %.2f, resting HR %.0f bpm, resp %.0f/min\n",
              x$fitness, x$effort, x$movement$strength, x$resting_hr, x$resp_rate))
  invisible(x)
}

#' @export
print.mm_intervention <- function(x, ...) {
  cat(sprintf("Multi-modal intervention (+%.2f fitness):\n", x$fitness_gain))
  cat(sprintf("  ROM %.3f -> %.3f (%+.0f%%),  mean HR %.0f -> %.0f bpm (%+.0f)\n",
              x$before["rom"], x$after["rom"], 100 * x$change["rom"] / x$before["rom"],
              x$before["mean_hr"], x$after["mean_hr"], x$change["mean_hr"]))
  invisible(x)
}
