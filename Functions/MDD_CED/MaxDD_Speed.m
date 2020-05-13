function [MaxDD, MaxDDIndex, MaxDDRecovery, RecoveryLength, DrawdownLength]...
    = MaxDD_Speed(Process,Path)

% Compute price proxy
CumuReturn = cumprod(1+Process);

% Computing dimension
[nDays,~] = size(Process);

% Computing Drawdowns
MaxDD = zeros(nDays-Path,1);
MaxDDIndex = zeros(nDays-Path,2);

for k = 1:nDays-Path
    
    [MaxDD(k), MaxDDIndex(k,:)] = maxdrawdown(CumuReturn(k:k+Path));
    MaxDDIndex(k,:) = MaxDDIndex(k,:) + k;
end

% Finding the recovery
MaxDDRecovery = zeros(nDays-Path,1);

for k = 1:size(MaxDD,1)
    
  EqualPeaks = find(CumuReturn(MaxDDIndex(k,2)+1:end) >= CumuReturn(MaxDDIndex(k,1)));
  
  if isempty(EqualPeaks)
      
      % If there is no recovery, index = NaN
      MaxDDRecovery(k) = NaN;
      
  else 
      
      % Find first recovery
      MaxDDRecovery(k) = EqualPeaks(1) + MaxDDIndex(k,2)+ 1;
       
  end
    
end

% Computing Speed
RecoveryLength = MaxDDRecovery - MaxDDIndex(:,2);
DrawdownLength = MaxDDIndex(:,2) - MaxDDIndex(:,1);

end