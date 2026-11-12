% Reproduces the analyses, tables and figures of the paper. Output goes to results/.
setup
cfg = paper_config();
out = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~isfolder(fullfile(out, 'figures')), mkdir(fullfile(out, 'figures')); end

%% Individual level (Figures 4 and 5)
results = run_individual_level(cfg);
results.rates = rate_statistics(results.rates, cfg.alpha);
save(fullfile(out, 'individual_level.mat'), '-struct', 'results');
writetable(results.rates, fullfile(out, 'individual_level_rates.csv'));

plot_rate_bars(results.rates, {'ttest', 'ensemble', 'ols'}, ...
    {'Two sample t-test/two group ANOVA', 'Ensemble t-test', 'OLS GLM'}, 0.05, ...
    fullfile(out, 'figures', 'Figure_4.png'));
plot_rate_bars(results.rates, {'arirls', 'bootstrap_ttest', 'bootstrap_ols'}, ...
    {'AR-IRLS GLM', 'Bootstrap t-test', 'Bootstrap OLS GLM'}, 0.5, ...
    fullfile(out, 'figures', 'Figure_5.png'));

%% Group level (Figure 6)
group = run_group_level(cfg);
save(fullfile(out, 'group_level.mat'), '-struct', 'group');
writetable(group.rates, fullfile(out, 'group_level_rates.csv'));

n_tests = cfg.group.n_replicates * cfg.group.channels_per_subj;
C = compare_group_rates(group.ols_fnr.rates, group.arirls_fnr.rates, n_tests, group.n_subjects, cfg.alpha);
writetable(C, fullfile(out, 'group_comparison_FNR.csv'));
disp(C)
plot_group_rates(group, 'FNR', fullfile(out, 'figures', 'Figure_6.png'));
