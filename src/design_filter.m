function [b, a] = design_filter(name, cfg)
% Butterworth coefficients for 'unfiltered' (empty), 'lowpass', 'highpass' or 'bandpass'.
nyq = cfg.fs / 2;
switch name
    case 'unfiltered'
        b = []; a = [];
    case 'lowpass'
        [b, a] = butter(cfg.filter_order, cfg.lowpass_cutoff / nyq, 'low');
    case 'highpass'
        [b, a] = butter(cfg.filter_order, cfg.highpass_cutoff / nyq, 'high');
    case 'bandpass'
        [b, a] = butter(cfg.filter_order, [cfg.highpass_cutoff cfg.lowpass_cutoff] / nyq);
end
end
