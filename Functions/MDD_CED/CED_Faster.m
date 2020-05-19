function [CondExp] = CED_Faster(Process, Path, Alpha)
%UNTITLED2 Summary of this function goes here

% Compute price proxy
CumuReturn = cumprod(1+Process);

% Computing dimension
[nDays,~] = size(Process);

% Computing Drawdowns
MaxDD = zeros(nDays-Path,1);

for k = 1:nDays-Path
    [MaxDD(k), ~] = maxdrawdown(CumuReturn(k:k+Path));
end

% Computing the conditional expected Drawdown
CondExp = mean(MaxDD(MaxDD>=quantile(MaxDD,1-Alpha)));

end

