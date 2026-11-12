function fig = plot_group_rates(group, measure, outfile)
% Group-level FNR (or FPR) vs number of subjects for OLS and AR-IRLS, Clopper-Pearson 95% CIs.
if nargin < 2, measure = 'FNR'; end
arms = {['ols_' lower(measure)], ['arirls_' lower(measure)]};
labels = {'OLS (bandpass filtered)', 'AR-IRLS (unfiltered)'};
markers = {'o-', 's-'};
cols = [0.03 0.19 0.42; 0.85 0.33 0.10];
n_subj = group.n_subjects;
n_tests = group.cfg.group.n_replicates * group.cfg.group.channels_per_subj;

fig = figure('Units', 'inches', 'Position', [1 1 3.3 2.6], 'Color', 'w');
ax = axes(fig); hold(ax, 'on');
for a = 1:2
    m = mean(group.(arms{a}).rates, 2);
    k = round(m * n_tests);
    [lo, hi] = clopper_pearson(k, n_tests * ones(size(k)));
    errorbar(ax, n_subj, m, m - lo(:), hi(:) - m, markers{a}, 'Color', cols(a, :), ...
        'MarkerFaceColor', cols(a, :), 'MarkerSize', 4, 'LineWidth', 1.2, 'CapSize', 3, 'DisplayName', labels{a});
end
if strcmp(measure, 'FPR')
    yline(ax, 0.05, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1, 'HandleVisibility', 'off');
    ylabel(ax, 'False Positive Rate');
else
    ylabel(ax, 'False Negative Rate');
end
xlabel(ax, 'Number of Subjects');
xlim(ax, [min(n_subj) - 0.5, max(n_subj) + 0.5]);
ax.FontSize = 8.5; ax.Box = 'off'; ax.LineWidth = 0.6;
legend(ax, 'Location', 'northeast', 'FontSize', 8, 'Box', 'off');
hold(ax, 'off');
if nargin >= 3 && ~isempty(outfile)
    exportgraphics(fig, outfile, 'Resolution', 600, 'BackgroundColor', 'white');
end
end
