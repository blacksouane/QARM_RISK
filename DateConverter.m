function [Freq] = DateConverter(Date,Frequency)
%Daily date to weekly or monthly ones

% INPUT:

% 1. Process : Daily Returns
% 2. Frequency : A value between Daily and Weekly

% Output: 

% 1. Freq : A vector of the changed frequency process


if Frequency == 1
    
   Freq=zeros(int16(round(length(Date))/5),1);
   Transition = datenum(Date);
   count = 1;
   
   for i=1:5:length(Date)
       
       if length(Date) - i > 5
           
        Freq(count) = Transition(i); 
        count = count + 1;
        
       end
   end
   
   Freq = datetime(Freq,'ConvertFrom','datenum','InputFormat','dd-MMM-yyyy');
                
                
elseif Frequency == 2

   Freq=zeros(int16(round(length(Date))/20),1);
   Transition = datenum(Date);
   count = 1;
   
   
   for i=1:20:length(Date)
       
        if length(Date) - i > 20 
            
        Freq(count) = Transition(i); 
        count = count + 1;
        
        end
   end
   
   Freq = datetime(Freq,'ConvertFrom','datenum','InputFormat','dd-MMM-yyyy'); 
    
else 
   disp('Not a viable frequency')




end

