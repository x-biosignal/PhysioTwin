---
title: 'PhysioTwin: Mechanistic simulation, sensor emulation and data assimilation for physiological and movement digital twins in R'
tags:
  - R
  - digital twin
  - physiological signals
  - data assimilation
  - Bayesian inference
  - in-silico trial
authors:
  - name: Yusuke Matsui
    affiliation: 1
affiliations:
  - name: Department of Rehabilitation Medicine, Nagoya University, Japan
    index: 1
date: 18 August 2026
bibliography: paper.bib
---

# Summary

`PhysioTwin` is an R package that builds *digital twins* of physiological and
movement systems: subject-specific mechanistic models that are personalised from
real measurements and then used for in-silico experiments. It couples three layers
that are usually found only in separate, domain-specific tools. First, mechanistic
generators produce signals from physical and physiological first principles: a
single-joint limb driven by forward dynamics, a Jansen-Rit neural-mass
electroencephalogram [@JansenRit1995], a McSharry limit-cycle electrocardiogram
[@McSharry2003], a Windkessel arterial-pressure model with a closed-loop baroreflex,
a respiration and respiratory-sinus-arrhythmia generator, and a closed-loop
cardiovascular-respiratory model whose Mayer waves emerge from the delayed baroreflex
[@deBoer1987]. Second, a forward sensor-emulation layer turns those clean signals
into realistic recordings (inertial units, optical markers with occlusion, surface
electromyography with crosstalk, multi-lead electrocardiograms, analogue-to-digital
conversion). Third, a data-assimilation layer personalises the twins to recorded
data with a family of state-space estimators (unscented, extended, particle and
ensemble Kalman filters), likelihood-free and full Bayesian inference (rejection and
sequential Monte-Carlo approximate Bayesian computation [@Toni2009], Markov-chain
Monte-Carlo, particle MCMC and particle Gibbs [@Andrieu2010; @Lindsten2014]), a
Gaussian-process surrogate and Bayesian-optimisation calibration [@Jones1998], and
Fisher-information optimal experimental design.

On top of these, `PhysioTwin` provides applications that make a twin useful:
single- and population-level in-silico intervention (predicting the effect of a
change, with effect size, responder rate and statistical power), a clinical decision
loop that gates a twin on its validation before it may advise and scores candidate
interventions with propagated uncertainty, labelled synthetic-training-data
generation, and a verification-and-validation harness. The package is written in
base R (its only hard dependency is `stats`), is extensively unit-tested, and ships
reproducible real-data validation cases (an electrocardiogram model against MIT-BIH,
a neural-mass model against a resting electroencephalogram, and the closed-loop
cardiovascular model against Fantasia heart-rate variability [@Goldberger2000]).

# Statement of need

Mechanistic simulators of individual physiological subsystems are well established
--- for example ECGSYN for the electrocardiogram [@McSharry2003] or neural-mass
models for the electroencephalogram [@JansenRit1995] --- and R has mature tools for
*analysing* recorded signals, such as `RHRV` for heart-rate variability. What is
missing is an integrated, open framework that (i) simulates *multiple* physiological
and movement systems from physical principles, (ii) emulates the sensors that would
record them, and (iii) closes the loop by assimilating real data to obtain a
*personalised* model that can then be experimented on in silico. This
predict-personalise-experiment loop is the essence of a digital twin, and it is
increasingly requested in rehabilitation, sports and neurophysiology, where one wants
to ask "how would this patient's movement or physiological response change under this
intervention?" before delivering it.

`PhysioTwin` fills this gap in R. Crucially, it treats the *validity* of a twin as a
first-class concern: it distinguishes verifying an estimator on synthetic data from
validating a model against reality, gates in-silico predictions on a per-subject
goodness-of-fit and identifiability check, propagates parameter uncertainty into
every intervention score, and reconciles mechanistic predictions with external
evidence rather than overriding it. This honest-scope discipline, together with the
breadth of mechanistic generators and the modern assimilation toolkit (particle MCMC,
particle Gibbs, sequential Monte-Carlo ABC and Bayesian-optimisation calibration) in
a dependency-light package, makes `PhysioTwin` a practical foundation for in-silico
intervention design and physically-consistent synthetic-data generation. It is part
of a broader Bioconductor-compatible ecosystem for physiological signals but stands
alone: it has no dependency on the other packages.

# Acknowledgements

We thank the PhysioNet project for the openly available recordings used in the
validation cases.

# References
