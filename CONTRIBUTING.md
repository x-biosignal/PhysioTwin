# Contributing to PhysioTwin

Contributions, bug reports and questions are welcome.

## Reporting bugs and requesting features

Please open an issue at
<https://github.com/x-biosignal/PhysioTwin/issues>. For a bug, include a minimal
reproducible example, the output of `sessionInfo()`, and what you expected to happen.

## Contributing code

1. Fork the repository and create a branch for your change.
2. Follow the existing style: base R only (the sole hard dependency is `stats`),
   roxygen2 documentation with runnable `@examples`, and a `testthat` test for every
   new exported function.
3. Run `devtools::document()`, `devtools::test()` and `R CMD check` locally; they
   should pass cleanly. Long-running Monte-Carlo tests use `skip_on_cran()`.
4. Open a pull request describing the change and how you verified it.

## Scope and validity

`PhysioTwin` is research tooling, not a medical device. New mechanistic models and
estimators should be validated: verify an estimator on synthetic data, and, where
possible, validate a model against real data (see `inst/validation/` for the
reproducible-case format). Please keep the honest-scope discipline of the existing
code --- state clearly what a method does and does not establish.

## Code of conduct

By contributing you agree to abide by the
[Contributor Covenant](https://www.contributor-covenant.org/) code of conduct.
