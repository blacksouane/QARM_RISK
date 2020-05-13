function [Freq] = FrequencyConverter(Process, Frequency)
%This function changes the frequency of observation of the process

% INPUT:

% 1. Process : Daily Returns
% 2. Frequency : A value between Daily and Weekly

% Output: 

% 1. Freq : A vector of the changed frequency process


if Frequency == 1
    
        % Compute distribution at weekly Freqeuncy
        Freq = zeros(int16(round(length(Process))/5), 1);
        count = 1;
        
        for i = 1:5:length(Process)
            
           if length(Process) - i > 5 
           temp = cumprod(1+Process(i:i+4)) -1;
           Freq(count) = temp(5);
           count = count + 1;
           end
           
        end
                
elseif Frequency == 2

        % Compute distribution at weekly Freqeuncy
        Freq = zeros(int16(round(length(Process))/20), 1);
        count = 1;
        
        for i = 1:20:length(Process)
           
           if length(Process) - i > 20 
           temp = cumprod(1+Process(i:i+19)) -1;
           Freq(count) = temp(5);
           count = count + 1;
           end
           
        end    
    
else 
   disp('Not a viable frequency')


end

