function [WeightsOpti] = RiskParityAllocation(Weights,Returns,Target,LengthPath, Alpha)
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
[~, asset] = size(Returns);

%Setting the size of the output. 
WeightsOpti = zeros(1500+2*LengthPath, asset);

%Setting the unused parameters of the optimisation (matlab obligates it)
A = []; %No linear constraint
b = []; %No linear constraint

%Disabling useless Warnings during optimisation
warning ( 'off' , 'MATLAB:nearlySingularMatrix')

% *************************************
%Loop optimizing the weights for each rebalancing
% *************************************
disp('Optimisation is starting !')
% Setting the objective function 
fun = @(x) -sum(log(x)); 

% Setting optimisation bounds on weights
lb = [];
ub = [];

% Setting linear constraint
Aeq = [];
beq = [];
        
% Setting options
options = optimoptions(@fmincon,'Algorithm','sqp',...
'MaxIterations',1000000,'ConstraintTolerance',1.0000e-10, ...
'OptimalityTolerance',1.0000e-10,'MaxFunctionEvaluations',...
1000000,'display','none');

for i = 2*LengthPath+1:2*LengthPath+1500
        

        % Optimizing the weights 
  WeightsOpti(i, :) = fmincon(fun,Weights(i, :),A,b,Aeq,beq,lb,ub,...
            @(x)CED_Constraint(x,Returns,Target,LengthPath,i, Alpha)...
            ,options);
     
        % Going for the next rebalancing 
        disp(i)
end
    
    disp('Optimisation is finished !')
   
    WeightsOpti = WeightsOpti(2*LengthPath+1:end,:);
end

