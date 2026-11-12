function S = make_surrogates(signal, theta)
% Phase-randomized surrogates (columns) with the amplitude spectrum of signal.
signal = signal(:).';
L = numel(signal);
if mod(L, 2) == 1, signal = [signal mean(signal)]; end
L2 = numel(signal);
Y = fft(signal, L2);
new_fft = zeros(size(theta, 1), L2);
new_fft(:, 2:(L2/2))     = Y(2:(L2/2)) .* exp(1i .* theta);
new_fft(:, 1)            = Y(1);
new_fft(:, L2/2 + 1)     = Y(L2/2 + 1);
new_fft(:, (L2/2)+2:end) = conj(fliplr(new_fft(:, 2:L2/2)));
S = ifft(new_fft, [], 2)';
S = S(1:L, :);
end
