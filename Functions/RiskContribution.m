function [MCR, CED_Portfolio] = RiskContribution(Data, Weights, Alpha, Window,h)

[~, NumAsset] = size(Data);
Signal = [1,-1];
MCR = zeros(length(Data) - 2*Window, NumAsset);

 for Asset = 1:NumAsset
      fprintf('Asset Number %s of the %s total assets',string(Asset),string(NumAsset))
      ChangedWeights = Weights;
      ChangedWeights(Asset) = ChangedWeights(Asset) + Signal(1)*h;
      Temp = PortfolioCED(ChangedWeights, Data, Alpha, Window);
      ChangedWeights(Asset) = ChangedWeights(Asset) + 2*Signal(2)*h;  
      Temp2 = PortfolioCED(ChangedWeights, Data, Alpha, Window);
      MCR(:,Asset) = (Temp - Temp2)./(2*h);
 end
 
 CED_Portfolio = PortfolioCED(Weights, Data, Alpha, Window);
 
 disp('Operation Completed')
end

