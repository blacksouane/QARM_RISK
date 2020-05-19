function [MCR] = MCR_VOL(x,Returns, LengthVol)
%This function compute the marginal risk contribution of each asset 

% INPUT : 

% 1. x : Vector/Matrix of weights 
% 2. Returns : Vector / Matrix of returns
% 3. LengthVol : A numeric value, the number of days to compute the
% volatility

% OUTPUT : 

% 1. MCR : Marginal contribution to risk of the strategy


% Pre-allocating the output
MCR = zeros(size(x,1)-LengthVol, size(x,2));

% Loop computing the contribution to risk for each days - the number of
% days to compute the initial covariance matrix.

for days = 1:length(x)-LengthVol
    
    Position = days+LengthVol;
    % 1. Computing the covariance
    CovMat = cov(Returns(days:days+LengthVol, :));
    
    % 2. Compute the marginal contribution to risk (in % of the total vol.)
    MCR(days,:) = (x(Position, :)'.*(CovMat*x(Position, :)')/...
        (x(Position, :)*CovMat*x(Position, :)'))'.*100;
    MCR(days, :) = MCR(days, :)/sum(MCR(days, :));
end

