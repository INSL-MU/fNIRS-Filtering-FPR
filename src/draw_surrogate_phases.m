function theta = draw_surrogate_phases(n_surrogates, n_samples)
% Uniform random phases for make_surrogates: n_surrogates x (L/2 - 1), L = even length.
L = n_samples + mod(n_samples, 2);
theta = unifrnd(0, 2*pi(), n_surrogates, (L/2) - 1);
end
