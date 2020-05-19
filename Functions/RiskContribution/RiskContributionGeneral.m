function [MCR] = RiskContributionGeneral(Data, Weights, Alpha, Window,h)

% This functions computes the marginal contribution to risk of each asset
% in a portfolio where the risk measure is the CED. 

% By definiton, the MCR of asset i is the partial derivative of the portfolio risk
% measure by the Weights of asset i. The closed form solution is relatively
% involved mathematically and therefore the implementation is done by
% numerical differentiation. 

% The drawback of this function is that the estimation is relatively noisy
% since it relies on empirical partial derivatives. 

% The numerical partial derivative is implemented using the "Symetric
% Derivative": 

% f'(w_i) = lim (h -> 0) (f(w+h) + f(x-h))/2h
% Where w +/- h means that the weights of asset i is changed by h and the
% other are the same.


% INPUT : 

% 1. Data : Daily Returns (Every frequency of observations is ok)
% 2. Weights : Initial Weights, i.e weights of the strategy.
% 3. Alpha : The level at which to compute the CED
% 4. Window : The path length and length of window to compute the CED
% 5. h : Numerical differentation value. Typically use 0.001

% OUTPUT : 

% 1. MCR : Marginal Contribution to risk where the risk is the CED for each
%          asset.
% 2. CED_Portfolio : The vector of CED at the level of the input 



[NumDays, NumAsset] = size(Data);
Signal = [1,-1];
MCR = zeros(length(Data) - 2*Window, NumAsset);

 for Asset = 1:NumAsset
     fprintf('Asset Number %s of the %s total assets \n',string(Asset),string(NumAsset))
     for Days = 1:NumDays-2*Window
      disp(Days)
      ChangedWeights = Weights(Days, :);
      ChangedWeights(Asset) = ChangedWeights(Asset) + Signal(1)*h;
      Temp = PortfolioCEDGeneral(ChangedWeights, Data(Days:Days+2*Window,:), Alpha, Window);
      ChangedWeights(Asset) = ChangedWeights(Asset) + 2*Signal(2)*h;  
      Temp2 = PortfolioCEDGeneral(ChangedWeights, Data(Days:Days+2*Window,:), Alpha, Window);
      MCR(Days,Asset) = (Temp - Temp2)./(2*h);
     end
 end
 
 disp('Operation Completed')
end

