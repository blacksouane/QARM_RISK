function [CondExp, MaxDD] = CED_Bitcoin(Process, Path, Alpha)
%UNTITLED2 Summary of this function goes here

% Computing dimension
[nDays,~] = size(Process);

% Computing Drawdowns
MaxDD = zeros(nDays-Path,1);

for k = 1:nDays-Path
    if ~mod(k,10000)
    disp(k);
    end
    [MaxDD(k), ~] = maxdrawdown(Process(k:k+Path));
end

% Computing the conditional expected Drawdown
CondExp = mean(MaxDD(MaxDD>=quantile(MaxDD,1-Alpha)));

end