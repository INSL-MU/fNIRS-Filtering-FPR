function [p, B, t] = fit_ols(X, Y)
% OLS GLM of each column of Y on X; p, B, t are regressors x columns.
[B, se] = lscov(X, Y);
t = B ./ se;
p = 2 * tcdf(-abs(t), size(X, 1) - size(X, 2));
end
