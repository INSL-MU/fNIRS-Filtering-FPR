function fig = plot_rate_bars(rates, methods, titles, fnr_ymax, outfile)
% FPR (left) and FNR (right) bar charts, one row per method, with Clopper-Pearson 95% CIs.
noises  = {'white', 'pink', 'rest'};
filters = {'unfiltered', 'bandpass', 'lowpass', 'highpass'};
cats    = {'White Noise', 'Pink Noise', 'Resting State\newlineData'};
legs    = {'Unfiltered', 'Bandpass', 'Lowpass', 'Highpass'};
cols    = [0.03 0.19 0.42; 0.85 0.33 0.10; 0.93 0.69 0.13; 0.82 0.82 0.82];

fig = figure('Units', 'inches', 'Position', [1 1 6.75 7.6], 'Color', 'w');
t = tiledlayout(fig, numel(methods), 2, 'TileSpacing', 'loose', 'Padding', 'compact');
for r = 1:numel(methods)
    for c = 1:2
        if c == 1, measure = 'FPR'; else, measure = 'FNR'; end
        [d, k, n] = rate_matrix(rates, methods{r}, measure, noises, filters);
        ax = nexttile(t);
        b = bar(ax, d, 'grouped');
        for j = 1:numel(filters)
            b(j).FaceColor = cols(j, :); b(j).EdgeColor = 'k'; b(j).LineWidth = 0.5;
        end
        hold(ax, 'on');
        if c == 1, yline(ax, 0.05, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1); end
        [lo, hi] = clopper_pearson(k, n);
        for j = 1:numel(filters)
            errorbar(ax, b(j).XEndPoints, d(:, j), d(:, j) - lo(:, j), hi(:, j) - d(:, j), ...
                'k', 'LineStyle', 'none', 'LineWidth', 0.75, 'CapSize', 3);
        end
        hold(ax, 'off');
        ax.XTickLabel = cats; ax.XTickLabelRotation = 0; ax.TickLabelInterpreter = 'tex';
        ax.FontSize = 8.5; ax.Box = 'off'; ax.LineWidth = 0.6;
        title(ax, titles{r}, 'FontSize', 9.5, 'FontWeight', 'bold');
        if c == 1
            ylabel(ax, 'False Positive Rate', 'FontSize', 9); ylim(ax, [0 1]);
        else
            ylabel(ax, 'False Negative Rate', 'FontSize', 9); ylim(ax, [0 fnr_ymax]);
        end
        if r == numel(methods) && c == 2
            legend(ax, b, legs, 'Location', 'northeast', 'FontSize', 8, 'Box', 'off');
        end
    end
end
annotation(fig, 'textbox', [0.005 0.955 0.05 0.04], 'String', '(a)', 'FontSize', 11, 'FontWeight', 'bold', 'EdgeColor', 'none', 'Margin', 0);
annotation(fig, 'textbox', [0.505 0.955 0.05 0.04], 'String', '(b)', 'FontSize', 11, 'FontWeight', 'bold', 'EdgeColor', 'none', 'Margin', 0);
if nargin >= 5 && ~isempty(outfile)
    exportgraphics(fig, outfile, 'Resolution', 600, 'BackgroundColor', 'white');
end
end

function [d, k, n] = rate_matrix(rates, method, measure, noises, filters)
d = zeros(numel(noises), numel(filters)); k = d; n = d;
for i = 1:numel(noises)
    for j = 1:numel(filters)
        row = rates.method == method & rates.measure == measure & rates.noise == noises{i} & rates.filter == filters{j};
        d(i, j) = rates.rate(row); k(i, j) = rates.k(row); n(i, j) = rates.n(row);
    end
end
end
