function [output] = SimulatedProcess(NumSim,AutoCorrParameter)
%Simulated an AR1 Process for the Drawdown process

% Pre-allocating the output
output = zeros(NumSim,1);

% Simulating the error term
Eps = normrnd(0,0.01,[NumSim,1]);

% Computing the process
for i = 1:NumSim
    
    if i == 1
      
        output(i) = Eps(i);
        
    else
        
        output(i) = output(i-1)*AutoCorrParameter + Eps(i);
    
    end
end



end

