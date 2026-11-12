function [lo, hi] = clopper_pearson(k, n, alpha)
% Exact binomial 1-alpha confidence interval for k successes in n trials (arrays allowed).
if nargin < 3, alpha = 0.05; end
[~, ci] = binofit(k(:), n(:), alpha);
lo = reshape(ci(:, 1), size(k));
hi = reshape(ci(:, 2), size(k));
end
