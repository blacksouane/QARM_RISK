function [Quantile,CondExp] = CED(Distribution,Alpha)
%CED Computes the conditional expected Drawdown given : 

% INPUT : 
 % 1. The distribution of maximum drawdowns over a given path
 % 2. The level of alpha, a numeric value between 0 and 1. 
 
% Output : 

% 1. Alphas quantile of the distribution
% 2. Conditional expected Drawdown 


% Computing the quantile
Quantile = quantile(Distribution,1-Alpha);


% Computing the conditional expected Drawdown
CondExp = mean(Distribution(Distribution>Quantile));

if isnan(CondExp) % if all DD are the same during the period
    CondExp = Quantile;
end
        
        
end

