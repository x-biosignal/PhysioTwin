# PhysioTwin

**Digital-twin simulation and data assimilation for physiological and movement
systems.** A digital twin couples a *mechanistic* model (whose parameters are
physical/clinical quantities) with *realistic sensor emulation* and *data
assimilation*, so a subject-specific model can be personalised from real
measurements and then used for in-silico experiments — predicting how movement or
physiology would change under an intervention, and generating richly-labelled
synthetic training data.

> **Status: prototype.** This release ships the **movement vertical slice** — a
> full end-to-end pipeline for one modality — proving the four-layer architecture.
> Other modalities (cardiac, neural, respiratory, network) are on the roadmap
> below. It is a private prototype (not published to the public org yet).

## The four layers

```
   parameters (strength, damping, stiffness, ...)
        │  1. MECHANISTIC SIMULATION   limbTwin() -> simulateTwin()
        ▼        (forward dynamics, RK4)
   true state (angle, velocity, acceleration)
        │  2. SENSOR EMULATION          imuMeasure()  (bias, noise, drift, scale)
        ▼
   realistic IMU stream  ◄──── real recorded data
        │  3. DATA ASSIMILATION         personalizeTwin()  (optim / UKF)
        ▼        + unscentedKalmanFilter()
   personalised twin  (parameters fit to the subject)
        │  4a. IN-SILICO INTERVENTION   insilicoIntervention()
        │  4b. TRAINING-DATA FACTORY    generateTrainingData()
        ▼
   predicted movement change  /  labelled synthetic dataset
```

## End-to-end example

```r
library(PhysioTwin)

# 1-2. a "real" subject, recorded through a noisy IMU
subject <- limbTwin(strength = 5.5, damping = 0.7, stiffness = 3.2)
rec <- imuMeasure(simulateTwin(subject, duration = 10, dt = 0.01, theta0 = 0.1),
                  imuSensor(gyro_bias = 0.015, gyro_noise = 0.035, accel_noise = 0.12))

# 3. personalise the twin from the IMU (recovers strength 5.5, damping 0.7, bias)
twin <- personalizeTwin(rec, limbTwin(strength = 4, damping = 0.3, stiffness = 3.2),
                        estimate = c("strength", "damping"), dt = 0.01)

# 4a. in-silico rehab: +40% strength, -30% damping -> predicted movement change
insilicoIntervention(twin$twin, scale = list(strength = 1.4, damping = 0.7))

# 4b. a labelled synthetic training set for ML (domain-randomised sensors)
generateTrainingData(200, seed = 7)
```

## Validity assessment (verification & validation)

Simulation validity is not the same as matching one dataset: it splits into
**verification** (are the equations solved correctly?) and **validation** (are
they the right equations for reality?), and a twin adds **assimilation
credibility** (does the personalised twin predict data it was not fitted on?).
`validateTwin()` runs the applicable checks and returns a structured report:

```r
validateTwin(limbTwin(strength = 6, damping = 0.5, stiffness = 6),
             estimate = c("strength", "damping"))
#> L0 verification : energy drift 9.6e-05, RK4 order 3.97      (numerics correct)
#> L1 identifiab.  : rank 2/2, condition 7.5                   (params recoverable)
#> L2 recovery     : strength bias +0.00 cover 100%; ...       (estimator unbiased)
#> L3 sensitivity  : strength=0.66 damping=6.32                (which params matter)
#> L4 prediction   : cross-condition rel-RMSE 0.001, face 2/2  (generalises + behaves)
#> L5 domain       : amplitude 3.46 rad (large-amplitude ...)  (operating envelope)
```

Validation against a subject's **real** data is done separately: personalise the
twin, then predict a held-out recording and quantify the error (the empirical /
cross-condition validation the harness scaffolds on synthetic data).

## Data assimilation

- `personalizeTwin(method = "optim")` — a batch/variational fit of the whole
  record (robust to strong non-linearity; the default), co-estimating a
  gyroscope-bias nuisance parameter.
- `personalizeTwin(method = "ukf")` — the recursive unscented Kalman filter
  (online; best in the mildly non-linear regime).
- `unscentedKalmanFilter()` — a stand-alone general nonlinear state estimator.

## Roadmap (comprehensive digital twin)

**Layer 1 — more mechanistic simulators** (reuse the ecosystem's Hill/EMG-driven
MSK, OpenSim, HD-EMG, EEG forward-head, EDA; add the missing generators): ECGSYN
cardiac ODE + IPFM tachogram + Windkessel/baroreflex, Jansen-Rit neural mass,
respiration + RSA, Kuramoto / central-pattern-generator coordination, network-
physiology coupling.

**Layer 2 — sensor emulation** beyond IMU: marker dropout/occlusion, EMG
electrode/amplifier (crosstalk, powerline, bandwidth, motion artifact), ECG lead,
PPG, wearable accelerometer; quantisation/sampling artefacts.

**Layer 3 — assimilation** beyond UKF/optim: EKF, particle filter / SMC, ensemble
Kalman filter, RTS smoother; Bayesian calibration (ABC / MCMC) with posteriors
over parameters; GP surrogate/emulator + global sensitivity analysis.

**Layer 4 — applications**: multi-parameter intervention design across modalities
(movement + physiology), population synthesis, and domain-randomised dataset
generation for AI.
