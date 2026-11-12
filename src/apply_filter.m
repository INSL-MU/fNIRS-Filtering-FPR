function y = apply_filter(x_padded, b, a, pad)
% Zero-phase filter of padded data, padding removed. Empty b: unfiltered.
if isempty(b)
    y = unpad(x_padded, pad);
else
    y = unpad(filtfilt(b, a, x_padded), pad);
end
end
