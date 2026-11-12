function cfg = paper_config()
% Parameters and random seeds of all simulations in the paper.

cfg.fs            = 3.125;   % Hz
cfg.n_samples     = 938;     % 300 s
cfg.n_channels    = 4000;    % records per data type
cfg.stim_period   = 20;      % s on, s off (7 blocks)
cfg.response_peak = 0.8;     % peak of the HRF response; noise has unit variance
cfg.pink_exponent = 0.8;     % pink noise P(f) ~ f^-0.8

cfg.filter_order    = 2;     % Butterworth, applied with filtfilt
cfg.lowpass_cutoff  = 0.1;   % Hz
cfg.highpass_cutoff = 0.01;  % Hz
cfg.filters         = {'unfiltered', 'bandpass', 'lowpass', 'highpass'};
cfg.pad_samples     = ceil(cfg.fs / cfg.highpass_cutoff);   % mirror padding for filtfilt

cfg.alpha        = 0.05;
cfg.ar_max_order = ceil(4 * cfg.fs);
cfg.n_bootstrap  = 1000;

cfg.group.n_subjects        = 2:10;
cfg.group.n_replicates      = 100;
cfg.group.channels_per_subj = 100;
cfg.group.response_rms      = 0.2;   % group simulations: unit-variance response times this

cfg.seed.data             = 1;
cfg.seed.bootstrap_ols    = 4;
cfg.seed.bootstrap_ttest  = 7;
cfg.seed.group_ols_fpr    = 100;
cfg.seed.group_ols_fnr    = 101;
cfg.seed.group_arirls_fpr = 102;
cfg.seed.group_arirls_fnr = 103;

cfg.resting_state_file = getenv('FNIRS_FILTERING_FPR_RESTDATA');
if isempty(cfg.resting_state_file)
    cfg.resting_state_file = fullfile(fileparts(mfilename('fullpath')), 'RestData.mat');
end
end
