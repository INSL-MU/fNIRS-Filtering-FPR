function results = run_individual_level(cfg, methods)
% FPR and FNR of each method on every data type and filter condition.
% results.rates: table (method, measure, noise, filter, k, n, rate); results.H: per-record decisions.
if nargin < 2
    methods = {'ttest', 'ensemble', 'ols', 'arirls', 'bootstrap_ttest', 'bootstrap_ols'};
end
[task, regressor, response] = make_task_design(cfg);
data = generate_datasets(cfg, response);
noises = {'white', 'pink', 'rest'};
X = [ones(cfg.n_samples, 1), (regressor - mean(regressor)) ./ std(regressor)];

H = struct();
for m = intersect({'ttest', 'ensemble', 'ols', 'arirls'}, methods, 'stable')
    fprintf('%s ...\n', m{1}); t0 = tic;
    for f = cfg.filters
        for d = noises
            for kind = {'null', 'task'}
                Y = data.(f{1}).(d{1}).(kind{1});
                switch m{1}
                    case 'ttest',    h = ttest_on_off(Y, task, cfg.alpha);
                    case 'ensemble', h = ensemble_ttest(Y, task, cfg.fs, cfg.alpha);
                    case 'ols',      p = fit_ols(X, Y); h = p(2, :) < cfg.alpha;
                    case 'arirls',   h = ar_irls_all(X, Y, cfg.ar_max_order, cfg.alpha);
                end
                H.(m{1}).(f{1}).(d{1}).(kind{1}) = logical(h(:)');
            end
        end
    end
    fprintf('  %.0f s\n', toc(t0));
end
if any(strcmp(methods, 'bootstrap_ols'))
    fprintf('bootstrap OLS ...\n'); t0 = tic;
    H.bootstrap_ols = bootstrap_ols_all(X, data, cfg);
    fprintf('  %.0f s\n', toc(t0));
end
if any(strcmp(methods, 'bootstrap_ttest'))
    fprintf('bootstrap t-test ...\n'); t0 = tic;
    H.bootstrap_ttest = bootstrap_ttest_all(data, task, cfg);
    fprintf('  %.0f s\n', toc(t0));
end
results.rates = summarize_rates(H, cfg.n_channels);
results.H = H;
results.cfg = cfg;
end

function data = generate_datasets(cfg, response)
% data.(filter).(noise).null / .task, padding removed. rest_diff holds the
% first-differenced resting-state records used by the bootstrap tests.
pad = cfg.pad_samples;
n_total = cfg.n_samples + 2*pad;
rng(cfg.seed.data);
padded.white = randn(n_total, cfg.n_channels);
padded.pink  = make_pink_noise(n_total, cfg.n_channels, cfg.fs, cfg.pink_exponent, pad);
padded.rest  = mirror_pad(load_resting_state(cfg), pad);
padded_response = zeros(n_total, 1);
padded_response(pad+1:end-pad) = response;
data = struct();
for f = cfg.filters
    [b, a] = design_filter(f{1}, cfg);
    for d = {'white', 'pink', 'rest'}
        data.(f{1}).(d{1}).null = apply_filter(padded.(d{1}), b, a, pad);
        data.(f{1}).(d{1}).task = apply_filter(padded.(d{1}) + padded_response, b, a, pad);
    end
    for kind = {'null', 'task'}
        x = padded.rest;
        if strcmp(kind{1}, 'task'), x = x + padded_response; end
        if ~isempty(b), x = filtfilt(b, a, x); end
        data.(f{1}).rest_diff.(kind{1}) = unpad([zeros(1, size(x, 2)); diff(x)], pad);
    end
end
end

function h = ar_irls_all(X, Y, max_order, alpha)
p = zeros(1, size(Y, 2));
parfor i = 1:size(Y, 2)
    pp = fit_ar_irls(X, Y(:, i), max_order);
    p(i) = pp(2);
end
h = p < alpha;
end

function H = bootstrap_ols_all(X, data, cfg)
% Surrogate phases are drawn once per record and test, in the order of the
% original runs (needed for seeded results to match the paper); each
% "filtered" draw is applied to the bandpass, lowpass and highpass versions.
% Resting-state records enter first-differenced, with a differenced regressor.
slots = {
    'white', 'unfiltered', 'null';  'pink', 'unfiltered', 'null';
    'white', 'filtered',   'null';  'pink', 'filtered',   'null';
    'white', 'unfiltered', 'task';  'pink', 'unfiltered', 'task';
    'white', 'filtered',   'task';  'pink', 'filtered',   'task';
    'rest',  'unfiltered', 'null';  'rest', 'filtered',   'null';
    'rest',  'unfiltered', 'task';  'rest', 'filtered',   'task'};
filtered = setdiff(cfg.filters, {'unfiltered'}, 'stable');
XTXi = pinv(X'*X);
X_diff = X; X_diff(:, 2) = [0; diff(X(:, 2))];
XTXi_diff = pinv(X_diff'*X_diff);
H = empty_H(cfg);
rng(cfg.seed.bootstrap_ols);
for i = 1:cfg.n_channels
    for s = 1:size(slots, 1)
        [noise, which, kind] = slots{s, :};
        theta = draw_surrogate_phases(cfg.n_bootstrap, cfg.n_samples);
        if strcmp(which, 'unfiltered'), flist = {'unfiltered'}; else, flist = filtered; end
        if strcmp(noise, 'rest')
            src = 'rest_diff'; Xs = X_diff; XTXis = XTXi_diff;
        else
            src = noise; Xs = X; XTXis = XTXi;
        end
        for f = flist
            y      = data.(f{1}).(src).(kind)(:, i);
            null_y = data.(f{1}).(src).null(:, i);
            H.(f{1}).(noise).(kind)(i) = bootstrap_ols(Xs, y, null_y, theta, cfg.alpha, XTXis);
        end
    end
    if mod(i, 500) == 0, fprintf('  %d / %d\n', i, cfg.n_channels); end
end
end

function H = bootstrap_ttest_all(data, task, cfg)
% One condition at a time over all records, in the order of the original runs.
% The two 'skip' conditions were run but not reported there; their draws are
% kept so that the later conditions line up with the paper.
% Resting-state records enter first-differenced.
conds = {
    'white', 'unfiltered', 'null';  'white', 'filtered', 'null';
    'pink',  'unfiltered', 'null';  'pink',  'filtered', 'null';
    'rest',  'unfiltered', 'null';  'rest',  'filtered', 'null';
    'rest',  'unfiltered', 'skip';  'rest',  'filtered', 'skip';
    'white', 'unfiltered', 'task';  'white', 'filtered', 'task';
    'pink',  'unfiltered', 'task';  'pink',  'filtered', 'task';
    'rest',  'unfiltered', 'task';  'rest',  'filtered', 'task'};
filtered = setdiff(cfg.filters, {'unfiltered'}, 'stable');
H = empty_H(cfg);
rng(cfg.seed.bootstrap_ttest);
for c = 1:size(conds, 1)
    [noise, which, kind] = conds{c, :};
    if strcmp(which, 'unfiltered'), flist = {'unfiltered'}; else, flist = filtered; end
    src = noise;
    if strcmp(noise, 'rest'), src = 'rest_diff'; end
    for i = 1:cfg.n_channels
        theta = draw_surrogate_phases(cfg.n_bootstrap, cfg.n_samples);
        if strcmp(kind, 'skip'), continue; end
        for f = flist
            y      = data.(f{1}).(src).(kind)(:, i);
            null_y = data.(f{1}).(src).null(:, i);
            H.(f{1}).(noise).(kind)(i) = bootstrap_ttest(y, task, null_y, theta, cfg.alpha);
        end
    end
    fprintf('  %s %s %s\n', which, noise, kind);
end
end

function H = empty_H(cfg)
H = struct();
for f = cfg.filters
    for d = {'white', 'pink', 'rest'}
        H.(f{1}).(d{1}).null = false(1, cfg.n_channels);
        H.(f{1}).(d{1}).task = false(1, cfg.n_channels);
    end
end
end

function T = summarize_rates(H, n)
rows = {};
for m = fieldnames(H)'
    for f = fieldnames(H.(m{1}))'
        for d = fieldnames(H.(m{1}).(f{1}))'
            h = H.(m{1}).(f{1}).(d{1});
            rows(end+1, :) = {m{1}, 'FPR', d{1}, f{1}, sum(h.null), n, sum(h.null) / n};   %#ok<AGROW>
            rows(end+1, :) = {m{1}, 'FNR', d{1}, f{1}, sum(~h.task), n, sum(~h.task) / n}; %#ok<AGROW>
        end
    end
end
T = cell2table(rows, 'VariableNames', {'method', 'measure', 'noise', 'filter', 'k', 'n', 'rate'});
for v = {'method', 'measure', 'noise', 'filter'}
    T.(v{1}) = categorical(T.(v{1}));
end
end
