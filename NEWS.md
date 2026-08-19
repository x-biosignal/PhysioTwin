# PhysioTwin 0.1.1

* **EEG 1/f aperiodic background.** `aperiodicNoise()` generates 1/f^exponent
  background noise, and `jansenRit(aperiodic=, aperiodic_exp=)` mixes it into the
  neural-mass column, turning its over-concentrated single-rhythm spectrum into a
  realistic broadband one: the relative alpha power drops from 0.99 to about 0.66,
  matching a real resting EEG, while the alpha peak survives (the dominant-rhythm
  search now runs above 3 Hz, over the 1/f floor). The Jansen-Rit validation case
  fits the 1/f amplitude so the model's relative alpha matches the recording.
* **Broadband low-frequency HRV model.** `cardioStateSpace(lf_model = "ar1")` and
  `personalizeCardioWaveform(lf_model = "ar1")` add a first-order-autoregressive
  (broadband, peakless) low-frequency component, so a record whose low-frequency
  variability is not a sharp Mayer rhythm is fitted without positing an arbitrary
  resonance frequency.

# PhysioTwin 0.1.0

First release: a digital-twin engine coupling mechanistic simulation of movement and
physiology with realistic sensor emulation and data assimilation.

* **Mechanistic generators.** A parameterised single-joint limb (forward dynamics)
  and mechanistic physiological/network generators: a Jansen-Rit neural-mass EEG
  column (`jansenRit()`), a McSharry limit-cycle electrocardiogram (`ecgsyn()`), a
  Windkessel arterial pressure model with a closed-loop baroreflex (`windkessel()`,
  `baroreflex()`), a respiration and respiratory-sinus-arrhythmia generator
  (`respiration()`, `rsaTachogram()`), Kuramoto and Matsuoka rhythms (`kuramoto()`,
  `cpgMatsuoka()`), and a closed-loop cardiovascular-respiratory model with emergent
  Mayer waves (`cardioRespiratory()`).
* **Sensor emulation.** Forward inertial-measurement-unit, optical-marker,
  surface-electromyography and electrocardiogram models, with richer chains for
  multi-channel EMG crosstalk (`emgArray()`), multi-lead ECG (`ecgMeasureLeads()`)
  and analogue-to-digital conversion (`adcQuantize()`).
* **Data assimilation.** Unscented, extended, particle and ensemble Kalman filters;
  rejection and sequential Monte-Carlo approximate Bayesian computation
  (`abcCalibration()`, `abcSMC()`); Markov-chain Monte-Carlo (`metropolis()`);
  particle MCMC and particle Gibbs for whole-waveform fitting (`particleMCMC()`,
  `particleGibbs()`); a Gaussian-process surrogate (`gpEmulator()`) and
  GP-accelerated calibration (`gpCalibrate()`); and Fisher-information optimal
  experimental design (`optimalDesign()`).
* **Personalisation and applications.** Fit a movement twin to recorded IMU data
  (`personalizeTwin()`) or the closed-loop cardiovascular twin to heart-rate/blood-
  pressure variability (`personalizeCardio()`, `personalizeCardioWaveform()`);
  single- and population-level in-silico intervention (`insilicoIntervention()`,
  `populationIntervention()`); a clinical decision loop with a validation gate,
  uncertainty-propagated scoring and evidence reconciliation (`clinicalDecision()`);
  and labelled synthetic training-data generation (`generateTrainingData()`).
* **Validity.** A verification-and-validation harness (`validateTwin()`,
  `profileLikelihood()`, `sobolIndices()`) and real-data validation cases (ECGSYN vs
  MIT-BIH, Jansen-Rit vs a resting EEG, and the closed loop vs Fantasia heart-rate
  variability) under `inst/validation/`.
