function [WeightsOpti] = RiskParity(Weights,Returns,Target,LengthVol)
%Optimisation of the Risk parity Weighting Scheme

%   This function take the following Inputs:

% *************************************
% INPUT
% *************************************

% Weights : A matrix of Initial Weights of assets X time dimension 
% Returns : A matrix of returns of assets X time dimension 
% Target :  A volatility target to compute the volatility constraint. 

% This yield the following output : 

% *************************************
% OUTPUT
% *************************************

% WeightsOpti : A matrix of assets X number of position of optimal weights
%               through the risk parity weighting scheme.

% *************************************
% Function
% *************************************

%Setting some parameters allowing the reproductibility of the function[
[NumDays, asset] = size(Returns);

%Setting the size of the output. 
WeightsOpti = zeros(NumDays-LengthVol, asset);

%Setting the unused parameters of the optimisation (matlab obligates it)
A = []; %No linear constraint
b = []; %No linear constraint

% Setting optimisation bounds on weights
lb = [];
ub = [];

% Setting linear constraint
Aeq = [];
beq = [];
        
%Disabling useless Warnings during optimisation
warning ( 'off' , 'MATLAB:nearlySingularMatrix')

% Setting the objective function 
fun = @(x) -sum(log(x)); 

% Setting options
options = optimoptions(@fmincon,'Algorithm','sqp',...
'MaxIterations',1000000,'ConstraintTolerance',1.0000e-10, ...
'OptimalityTolerance',1.0000e-10,'MaxFunctionEvaluations',...
1000000,'display','none');
position = 1;


% *************************************
%Loop optimizing the weights for each rebalancing
% *************************************
disp('Optimisation is starting !')
for i = LengthVol+1:NumDays
        
    %Computing the covariance matrix
    CovMat = cov(Returns(i-LengthVol:i,:));
    
    % Optimizing the weights 
    x = fmincon(fun,Weights(i, :)...
    ,A,b,Aeq,beq,lb,ub,...
    @(x)VolConstraint(x,Target,CovMat),options);
     
    % Rescaling the weights
    WeightsOpti(position, :) = x/sum(x);
    
    % Going for the next rebalancing 
     position = position + 1;
     disp(position)
end
    
    disp('Optimisation is finished !')

end

