function H = ttest_on_off(Y, task, alpha)
% Welch t-test of task-on vs task-off samples, per column of Y.
on = task(:, 1) == 1;
H = logical(ttest2(Y(on, :), Y(~on, :), 'Alpha', alpha, 'Vartype', 'unequal'));
end
