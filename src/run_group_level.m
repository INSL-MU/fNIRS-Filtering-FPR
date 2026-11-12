function group = run_group_level(cfg, arms)
% Group-level (random effects) FPR and FNR for 2..10 subjects of pink-noise channels.
% Arms: ols_fpr, ols_fnr (bandpass + OLS), arirls_fpr, arirls_fnr (unfiltered + AR-IRLS).
% group.(arm).rates is group sizes x replicates; group.rates a summary table.
if nargin < 2, arms = {'ols_fpr', 'ols_fnr', 'arirls_fpr', 'arirls_fnr'}; end
spec.ols_fpr    = struct('method', 'ols',    'filter', 'bandpass',   'task', false, 'seed', cfg.seed.group_ols_fpr);
spec.ols_fnr    = struct('method', 'ols',    'filter', 'bandpass',   'task', true,  'seed', cfg.seed.group_ols_fnr);
spec.arirls_fpr = struct('method', 'arirls', 'filter', 'unfiltered', 'task', false, 'seed', cfg.seed.group_arirls_fpr);
spec.arirls_fnr = struct('method', 'arirls', 'filter', 'unfiltered', 'task', true,  'seed', cfg.seed.group_arirls_fnr);

[~, regressor] = make_task_design(cfg);
X = [ones(cfg.n_samples, 1), (regressor - mean(regressor)) ./ std(regressor)];
pad = cfg.pad_samples;
n_total = cfg.n_samples + 2*pad;
padded_response = zeros(n_total, 1);
padded_response(pad+1:end-pad) = cfg.group.response_rms * regressor;
n_ch = cfg.group.channels_per_subj;
sizes = cfg.group.n_subjects;

group = struct('n_subjects', sizes, 'cfg', cfg);
rows = {};
for a = arms
    sp = spec.(a{1});
    [b, aa] = design_filter(sp.filter, cfg);
    fprintf('%s ...\n', a{1}); t0 = tic;
    rates = zeros(numel(sizes), cfg.group.n_replicates);
    rng(sp.seed);
    for s = 1:numel(sizes)
        for r = 1:cfg.group.n_replicates
            Y = zeros(cfg.n_samples, n_ch, sizes(s));
            for subj = 1:sizes(s)                       % sequential: consumes the random stream
                pink = make_pink_noise(n_total, n_ch, cfg.fs, cfg.pink_exponent, pad);
                if sp.task, pink = pink + padded_response; end
                Y(:, :, subj) = apply_filter(pink, b, aa, pad);
            end
            betas = first_level_betas(X, Y, sp.method, cfg.ar_max_order);
            h = ttest(betas', 0, 'Alpha', cfg.alpha);    % one-sample t-test across subjects, per channel
            if sp.task, rates(s, r) = sum(~h) / n_ch; else, rates(s, r) = sum(h) / n_ch; end
        end
        fprintf('  n = %2d: %.4f  (%.0f s)\n', sizes(s), mean(rates(s, :)), toc(t0));
        if sp.task, measure = 'FNR'; else, measure = 'FPR'; end
        rows(end+1, :) = {a{1}, sp.method, measure, sizes(s), mean(rates(s, :)), std(rates(s, :))}; %#ok<AGROW>
    end
    group.(a{1}).rates = rates;
    group.(a{1}).spec = sp;
end
group.rates = cell2table(rows, 'VariableNames', {'arm', 'method', 'measure', 'n_subjects', 'mean_rate', 'sd_rate'});
end

function betas = first_level_betas(X, Y, method, max_order)
[~, n_ch, n_subj] = size(Y);
switch method
    case 'ols'
        betas = zeros(n_ch, n_subj);
        for subj = 1:n_subj
            [~, B] = fit_ols(X, Y(:, :, subj));
            betas(:, subj) = B(2, :)';
        end
    case 'arirls'
        Yflat = reshape(Y, size(Y, 1), n_ch * n_subj);
        bflat = zeros(1, n_ch * n_subj);
        parfor j = 1:n_ch * n_subj
            [~, B] = fit_ar_irls(X, Yflat(:, j), max_order);
            bflat(j) = B(2);
        end
        betas = reshape(bflat, n_ch, n_subj);
end
end
