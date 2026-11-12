function xp = mirror_pad(x, pad)
% Extends each column by mirroring its first and last pad samples.
xp = [flipud(x(1:pad, :)); x; flipud(x(end-pad+1:end, :))];
end
