function [task, regressor, response] = make_task_design(cfg)
% Block-design boxcar, its standardized HRF regressor, and the response added to the noise.
time = 0:(1/cfg.fs):((cfg.n_samples-1)/cfg.fs);
task = mod(floor(time(:) / cfg.stim_period), 2);
h = canonical_hrf(1 / (time(2) - time(1)));
convolved = filter(h, 1, task);
regressor = (convolved - mean(convolved)) ./ std(convolved);
response  = cfg.response_peak * (convolved ./ max(convolved));
end
