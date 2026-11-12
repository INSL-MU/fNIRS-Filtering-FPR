function h = canonical_hrf(fs)
% Canonical double-gamma HRF sampled at fs over 0-32 s, unit sum (SPM / AnalyzIR defaults).
t = (0:1/fs:32)';
h = 1^4 * t.^3 .* exp(-t) / gamma(4) - (1/6) * 1^16 * t.^15 .* exp(-t) / gamma(16);
h = h / sum(h);
end
