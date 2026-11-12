function T = compare_group_rates(rates_A, rates_B, n_tests, n_subjects, alpha)
% Compares two independent group-simulation arms (levels x replicates matrices of rates)
% at each group size: pooled rates, risk difference with Newcombe CI, Cohen's h,
% Fisher's exact test raw and Holm-corrected over group sizes.
if nargin < 5, alpha = 0.05; end
kA = round(mean(rates_A, 2) * n_tests);
kB = round(mean(rates_B, 2) * n_tests);
pA = kA / n_tests;
pB = kB / n_tests;
p_raw = zeros(numel(n_subjects), 1);
for i = 1:numel(n_subjects)
    p_raw(i) = fisher_2x2(kA(i), n_tests, kB(i), n_tests);
end
[rd_lo, rd_hi] = newcombe_ci(kA, n_tests, kB, n_tests, alpha);
h = 2*asin(sqrt(pA)) - 2*asin(sqrt(pB));
T = table(n_subjects(:), kA, pA, kB, pB, pA - pB, rd_lo, rd_hi, h, p_raw, holm(p_raw), ...
    'VariableNames', {'n_subjects', 'k_A', 'rate_A', 'k_B', 'rate_B', 'RD', 'RD_lo', 'RD_hi', 'cohens_h', 'p_raw', 'p_holm'});
end

function p = fisher_2x2(x1, n1, x2, n2)
if x1 == x2 || (x1 == 0 && x2 == 0) || (x1 == n1 && x2 == n2)
    p = 1;
    return
end
[~, p] = fishertest([x1, n1 - x1; x2, n2 - x2]);
end

function [lo, hi] = wilson_ci(x, n, alpha)
z = norminv(1 - alpha/2);
phat = x ./ n;
denom = 1 + z^2 ./ n;
center = (phat + z^2 ./ (2*n)) ./ denom;
halfw = (z ./ denom) .* sqrt(phat .* (1 - phat) ./ n + z^2 ./ (4 * n.^2));
lo = max(0, center - halfw);
hi = min(1, center + halfw);
end

function [lo, hi] = newcombe_ci(x1, n1, x2, n2, alpha)
% Newcombe (1998) hybrid score interval for a difference of proportions.
[l1, u1] = wilson_ci(x1, n1, alpha);
[l2, u2] = wilson_ci(x2, n2, alpha);
p1 = x1 ./ n1; p2 = x2 ./ n2; d = p1 - p2;
lo = d - sqrt((p1 - l1).^2 + (u2 - p2).^2);
hi = d + sqrt((u1 - p1).^2 + (p2 - l2).^2);
end

function p_adj = holm(p)
p = p(:); m = numel(p);
[ps, idx] = sort(p);
adj = min(cummax(ps .* (m:-1:1)'), 1);
p_adj = zeros(m, 1);
p_adj(idx) = adj;
end
