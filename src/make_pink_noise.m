function noise = make_pink_noise(n_samples, n_series, fs, exponent, pad)
% n_samples x n_series noise with amplitude spectrum f^-exponent and random phases.
% Columns are standardized over samples pad+1:end-pad. n_samples must be even.
if nargin < 5, pad = 0; end
freqs = 0:fs/n_samples:fs/2;
amplitude = arrayfun(@(f) f^(-exponent), freqs);
phase = unifrnd(0, 2*pi(), (n_samples/2)+1, n_series);
S = repmat(amplitude', [1 n_series]) .* exp(1i*phase);
S(1,:) = 0;
S2 = [S(1,:)' S(2:end,:)'/2 fliplr(conj(S(2:end,:)'))/2];
noise = ifft(S2, [], 2)';
noise = noise(2:end, :);
core = noise(pad+1:end-pad, :);
noise = (noise - mean(core, 1)) ./ std(core);
end
