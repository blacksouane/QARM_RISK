function [CED_PortfolioDistribution,ES_PortfolioDistribution] = CED_portfolio(Asset,RiskF,wA,Alpha,Path)
%Compute the maximum drawdown distribution of the Process Input over a path
% of length length(Path)

%   Input: 

% 1. Asset : Timeseries of returns
% 2. RiskF : Timeseries of risk free rate
% 3. wA = weight of market porfolio
% 4. The level of alpha, a numeric value between 0 and 1. 
% 5. Path : Integer Number of days (all period are working)

% Output: 

% 1. CED_PortfolioDistribution : CED for the portfolio

%******************************************************************

% Number of days in the process
[nDays,~] = size(Asset);

% Transform Returns into price
Portfolio = zeros(nDays,1);

for i=1:nDays
    Portfolio(i)=wA*Asset(i)+(1-wA)*RiskF(i);
end

[MaxDis] = MDD_Distribution(Portfolio,Path);

[~,CED_PortfolioDistribution] = CED(MaxDis,Alpha);

[~, ES_PortfolioDistribution, ~] = ES(Portfolio,Alpha);

        
end

