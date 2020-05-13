function [VaR, ES, LossDistribution] = ES(Process,Alpha)

% Pre - Processing
[nDays, nAsset] = size(Process);

% Transform Returns into price
Price = ones(nDays+1,nAsset);
Price(2:end,:) = cumprod(1+Process);

% Computing the Var
Distribution = sort((Price(1:end-1,:) - Price(2:end,:))./Price(1:end-1,:));
LossDistribution = abs(Distribution(Distribution < 0));
VaR = quantile(LossDistribution,Alpha);

% Computing the Expected shortfall
ES = mean(LossDistribution(LossDistribution >= VaR));
end