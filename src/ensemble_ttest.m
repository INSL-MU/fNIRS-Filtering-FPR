function H = ensemble_ttest(Y, task, fs, alpha)
% Welch t-test on the block average of each column (epochs -10 to +30 s around onsets),
% task-block samples vs the rest.
onsets = find(diff(task(:, 1)) == 1) + 1;
pre  = ceil(10 * fs);
post = floor(20 * fs) + floor(10 * fs);
avg = zeros(pre + post + 1, size(Y, 2));
for k = 1:numel(onsets)
    avg = avg + Y(onsets(k)-pre : onsets(k)+post, :);
end
avg = avg ./ numel(onsets);
block = pre+1 : pre+floor(20*fs)+1;
rest  = [1:pre, pre+floor(20*fs)+2 : pre+post+1];
H = logical(ttest2(avg(block, :), avg(rest, :), 'Alpha', alpha, 'Vartype', 'unequal'));
end
