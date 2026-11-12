function [H, B, CI] = bootstrap_ols(X, y, null_y, theta, alpha, XTXi)
% OLS task coefficient of y tested against the coefficients of phase-randomized
% surrogates of null_y (central 1-alpha interval).
if nargin < 6, XTXi = pinv(X'*X); end
null_B = XTXi(2, :) * (X' * make_surrogates(null_y, theta));
B = XTXi * (X' * y);
B = B(2);
CI = [prctile(null_B, 100*alpha/2) prctile(null_B, 100*(1-alpha/2))];
H = ~((B < CI(2)) & (B > CI(1)));
end
