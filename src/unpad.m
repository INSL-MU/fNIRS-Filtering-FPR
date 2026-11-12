function x = unpad(xp, pad)
x = xp(pad+1:end-pad, :);
end
