function [p, B, ar_order] = fit_ar_irls(X, y, max_order)
% AR-IRLS GLM (Barker et al. 2013): AR prewhitening of the residuals (order by BIC) and
% bisquare robust refit, iterated until the coefficients settle.
warn_state = warning('off', 'stats:statrobustfit:IterationLimit');
cleanup = onCleanup(@() warning(warn_state));

B  = (X'*X) \ X' * y;
B0 = 1e6 * ones(size(B));
iter = 0;
while norm(B - B0) / norm(B0) > 1e-2 && iter < 10
    B0 = B;
    a = ar_fit_bic(y - X*B, max_order);      % [c; a1; ...; ap]
    f = [1; -a(2:end)];
    Xf = filter_keep_offset(f, X);
    yf = filter_keep_offset(f, y);
    [B, S] = robustfit(Xf, yf, 'bisquare', 4.685, 'off');
    iter = iter + 1;
end
covB = pinv(Xf'*Xf) * S.robust_s^2;
t = B ./ sqrt(diag(covB));
p = 2 * (1 - tcdf(abs(t), length(y) - 1));
ar_order = numel(a) - 1;
end

function out = filter_keep_offset(f, y)
% FIR filter after removing the first sample of each column, offset restored.
y1 = y(1, :);
out = filter(f, 1, y - y1) + sum(f) * y1;
end

function coef = ar_fit_bic(y, max_order)
% Least-squares AR fit on the series and its reverse, order 0..max_order chosen by BIC.
n = length(y);
max_order = max(min(max_order, n - 1), 1);
X  = [ones(2*n, 1) [lag_matrix(y, max_order); lag_matrix(flipud(y), max_order)]];
yy = [y; flipud(y)];
[Q, R] = qr(X, 0);
invR = pinv(R);
m = 2*n;
loglik = nan(size(X, 2), 1);
for k = 1:numel(loglik)
    b = invR(1:k, 1:k) * Q(:, 1:k)' * yy;
    r = yy - X(:, 1:k) * b;
    loglik(k) = -m/2 * log(2*pi*mean(r.^2)) - m/2;
end
bic = -2*loglik + (1:numel(loglik))' .* log(m);
[~, best] = min(bic);
coef = invR(1:best, 1:best) * Q(:, 1:best)' * yy;
end

function L = lag_matrix(y, max_lag)
n = numel(y);
L = zeros(n, max_lag);
for k = 1:max_lag
    L(k+1:end, k) = y(1:n-k);
end
end
