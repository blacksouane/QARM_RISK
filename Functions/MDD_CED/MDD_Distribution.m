function [MaxDis] = MDD_Distribution(Process,Path)
%Compute the maximum drawdown distribution of the Process Input over a path
% of length length(Path)

%   Input: 

% 1. Process : Timeseries of returns
% 2. Path : Integer Number of days (all period are working)

% Output: 

% 1. Distribution : Distribution of Maximum Drawdown over the Path

%******************************************************************

% Number of days in the process
[nDays,nAsset] = size(Process);

% Transform Returns into price
Price = ones(nDays+1,nAsset);
Price(2:end,:) = cumprod(1+Process);

% Preallocate Variable
MaxDis = zeros(nDays+1-Path,nAsset);

% Computing the distribution
    for k = 1:nDays+1-Path
    MaxData = Price(k,:);
    MaxDD = zeros(1,nAsset);
        for i = 1:Path
		MaxData = max(MaxData, Price(k+i,:));
		DD = (MaxData - Price(i+k,:)) ./ MaxData;
            if any(DD > MaxDD)
			p = DD > MaxDD;
			MaxDD(p) = DD(p);
            end
        end
    MaxDis(k,:) = MaxDD;
    end


        
end

