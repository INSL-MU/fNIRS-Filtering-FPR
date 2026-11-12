function T = rate_statistics(rates, alpha, nominal)
% Adds Clopper-Pearson intervals and, for FPR rows, a one-tailed exact binomial test of
% FPR > nominal, Bonferroni-corrected over all rows of the table.
if nargin < 2, alpha = 0.05; end
if nargin < 3, nominal = 0.05; end
T = rates;
[T.ci_lo, T.ci_hi] = clopper_pearson(T.k, T.n, alpha);
T.p_binomial = nan(height(T), 1);
T.inflated = false(height(T), 1);
fpr = T.measure == 'FPR';
p = 1 - binocdf(T.k(fpr) - 1, T.n(fpr), nominal);
T.p_binomial(fpr) = p;
T.inflated(fpr) = p < alpha / height(T);
end
