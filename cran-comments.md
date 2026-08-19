## Submission

This is a new submission (PhysioTwin 0.1.1).

## Test environments

* local: Ubuntu, R (release), `R CMD check --as-cran`
* (to add before submission) win-builder (devel and release), macOS builder,
  R-hub linux/windows/macOS

## R CMD check results

0 errors | 0 warnings | 0 notes (local `--as-cran`).

## Notes for CRAN

* The package suggests two packages hosted on an r-universe repository
  (`PhysioRehab`, `PhysioAppKit`), declared via `Additional_repositories:
  https://x-biosignal.r-universe.dev`. They are used only by an optional bridge
  (`rehabEvaluate()`), which is guarded by `requireNamespace()`; all tests and
  examples that use them are skipped when they are not installed, so the package
  checks and works fully without them.
* Long-running Monte-Carlo tests (particle filters, MCMC, Bayesian optimisation and
  Monte-Carlo recovery studies) are marked `skip_on_cran()`; the fast unit tests that
  verify the core behaviour run on CRAN.
* All examples run in a few seconds; the heavier ones are wrapped in `\donttest{}`.
* The package writes no files outside `tempdir()` and makes no network calls.
