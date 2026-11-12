function rest = load_resting_state(cfg)
% First cfg.n_samples samples of the first cfg.n_channels resting-state channels, standardized.
if ~isfile(cfg.resting_state_file)
    error('Resting-state file not found: %s (see README)', cfg.resting_state_file);
end
S = load(cfg.resting_state_file, 'restData');
rest = S.restData(1:cfg.n_samples, 1:cfg.n_channels);
rest = (rest - mean(rest, 1)) ./ std(rest, 0, 1);
end
