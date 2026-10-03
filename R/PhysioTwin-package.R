#' PhysioTwin: Digital-Twin Simulation and Data Assimilation for Physiology and Movement
#'
#' PhysioTwin builds subject-specific *digital twins*: it couples a mechanistic
#' model (whose parameters are physical or clinical quantities) with realistic
#' sensor emulation and data assimilation, so a model can be personalised from
#' real measurements and then used for in-silico experiments -- predicting how
#' movement or physiology would change under an intervention, and generating
#' richly labelled synthetic training data. The engine is organised as four
#' layers. The **movement** modality is worked end-to-end as the reference
#' vertical slice (simulate, sense, personalise, intervene / generate); the
#' physiological generators and estimators below are mechanistic and statistical
#' building blocks for coupling the other modalities into the same loop.
#'
#' @details
#' **1. Mechanistic simulators (the model)**
#'
#' * Movement limb -- [limbTwin()], [simulateTwin()]
#' * Cardiac -- [ecgsyn()], [windkessel()], [baroreflex()], [cardioStateSpace()],
#'   [cardioRespiratory()]
#' * Neural (EEG) -- [jansenRit()], [aperiodicNoise()], [bandPower()]
#' * Respiration / HRV -- [respiration()], [rsaTachogram()]
#' * Coordination rhythms -- [kuramoto()], [cpgMatsuoka()]
#' * Multimodal coupling -- [simulateMultimodal()], [multimodalTwin()]
#'
#' **2. Sensor emulation (turn a clean signal into a realistic recording)**
#'
#' * Inertial (IMU) -- [imuSensor()], [imuMeasure()]
#' * Optical marker -- [markerSensor()], [markerMeasure()], [fillMarkerGaps()]
#' * Electromyography -- [emgElectrode()], [emgMeasure()], [emgArray()],
#'   [emgArrayMeasure()], [emgEnvelope()], [crosstalkMatrix()]
#' * Electrocardiogram -- [ecgLead()], [ecgLeadSet()], [ecgMeasure()],
#'   [ecgMeasureLeads()], [detectRpeaks()]
#' * Analogue-to-digital -- [adcSample()], [adcQuantize()], [quantizationNoise()]
#'
#' **3. Data assimilation (personalise the twin to a subject's data)**
#'
#' * Fit a twin -- [personalizeTwin()], [personalizeCardio()],
#'   [personalizeCardioWaveform()]
#' * State estimators -- [unscentedKalmanFilter()], [extendedKalmanFilter()],
#'   [particleFilter()], [ensembleKalmanFilter()]
#' * Likelihood-free / Bayesian calibration -- [abcCalibration()], [abcSMC()],
#'   [metropolis()], [particleMCMC()], [particleGibbs()]
#' * Surrogate, sensitivity and design -- [gpEmulator()], [gpCalibrate()],
#'   [sobolIndices()], [profileLikelihood()], [optimalDesign()]
#'
#' **4. Applications**
#'
#' * In-silico intervention -- [insilicoIntervention()],
#'   [insilicoInterventionMultimodal()], [intervention()],
#'   [defaultInterventions()], [evaluateIntervention()], [scoreInterventions()],
#'   [populationIntervention()]
#' * Synthetic training data -- [generateTrainingData()]
#' * Validity / verification harness -- [validateTwin()], [twinReadiness()]
#' * Clinical decision loop -- [clinicalDecision()], [reconcileEvidence()],
#'   [rehabEvaluate()]
#'
#' **Where to go next**
#'
#' Start with `vignette("PhysioTwin")`, which runs the four movement layers
#' end-to-end ([limbTwin()] -> [imuMeasure()] -> [personalizeTwin()] ->
#' [insilicoIntervention()] and [generateTrainingData()]). The clinical bridge
#' [rehabEvaluate()] holds a twin's predictions to account against a PhysioRehab
#' single-case evaluation, and ordinal outcomes can be scaled with PhysioAppKit.
#'
#' @references
#' Jansen BH, Rit VG (1995) Biol Cybern 73:357-366 (neural-mass EEG).
#' McSharry PE, Clifford GD, Tarassenko L, Smith LA (2003) IEEE Trans Biomed Eng
#' 50:289-294 (ECG). Julier SJ, Uhlmann JK (1997) Proc SPIE 3068:182-193 (UKF).
#'
#' @keywords internal
"_PACKAGE"
