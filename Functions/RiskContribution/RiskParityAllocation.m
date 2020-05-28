function [WeightsOpti, MCR] = RiskParityAllocation(Weights,Returns,Target,LengthPath, Alpha)
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
MCR = zeros(1500+2*LengthPath, asset);
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
options = optimoptions(@fmincon,'Algorithm','interior-point',...
'MaxIterations',1000000,'ConstraintTolerance',1.0000e-6, ...
'OptimalityTolerance',1.0000e-6,'MaxFunctionEvaluations',...
1000000);

for i = 2*LengthPath+1:2*LengthPath+1500

      % Optimizing the weights
       
  [x,~,~,~,lambda,~,~] = fmincon(fun,Weights(i, :),A,b,Aeq,beq,lb,ub,...
          @(x)CED_Constraint(x,Returns, Target, LengthPath, i, Alpha)...
            ,options);
        
        
     WeightsOpti(i, :) = x;
     
     for assets = 1:length(x)
     MCR(i, assets) = 1/(x(assets)*lambda.ineqnonlin);
     end
     disp(MCR)
     % Going for the next rebalancing 
     disp(i)
end
    
    disp('Optimisation is finished !')
   
    % We delete the part of the matrix that was here for the computations
    WeightsOpti = WeightsOpti(2*LengthPath+1:end,:);
    MCR = MCR(2*LengthPath+1:end,:);
end

