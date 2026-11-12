function [H, p] = bootstrap_ttest(y, task, null_y, theta, alpha)
% Welch t statistic (task-on vs task-off) of y tested against the statistics of
% phase-randomized surrogates of null_y (central 1-alpha interval).
on = task(:, 1) == 1;
S = make_surrogates(null_y, theta);
[~, ~, ~, st] = ttest2(S(on, :), S(~on, :), 'Vartype', 'unequal');
null_t = st.tstat;
[~, ~, ~, st] = ttest2(y(on), y(~on), 'Vartype', 'unequal');
t = st.tstat;
CI = prctile(null_t, [100*alpha/2, 100*(1-alpha/2)]);
H = CI(1) > t | CI(2) < t;
p = 1 - sum(abs(null_t) < abs(t)) / numel(null_t);
end
