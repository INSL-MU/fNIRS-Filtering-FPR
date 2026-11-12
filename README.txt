fNIRS-Filtering-FPR

Code for: Zinos A, Wagner JC, Whelan HT, Beardsley SA. Reliable individual-level
fNIRS requires statistical control of serial correlation rather than frequency
filtering. Neurophotonics (submitted).

Simulates false positive and false negative rates of six analysis methods
(t-test, ensemble-averaged t-test, OLS GLM, AR-IRLS GLM, phase-randomized
bootstrap t-test and GLM) on white noise, pink noise and resting-state fNIRS
recordings, unfiltered and after lowpass, highpass and bandpass filtering,
plus a group-level simulation comparing OLS and AR-IRLS.

Requirements
  MATLAB R2023a, Signal Processing Toolbox, Statistics and Machine Learning
  Toolbox. Parallel Computing Toolbox is optional (AR-IRLS uses parfor).

Data
  The resting-state recordings (restData, samples x 4000, from Wagner et al.
  2025, NeuroImage 321:121492) are available from the corresponding author on
  request. Put RestData.mat in the repository folder or set the environment
  variable FNIRS_FILTERING_FPR_RESTDATA to its path.

Run
  setup
  run_paper

  Writes the rates (CSV and MAT) and Figures 4 to 6 to results/. Parameters
  and random seeds are in paper_config.m.
