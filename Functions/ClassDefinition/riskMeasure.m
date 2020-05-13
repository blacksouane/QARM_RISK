classdef riskMeasure
   
    properties 
        
       ConfidenceLevel {mustBeGreaterThanOrEqual(ConfidenceLevel,0),...
           mustBeLessThanOrEqual(ConfidenceLevel,1)} = 0.05
       PathLength = 125 % Setting default value
       MDD_Distribution
       MDD
       CED
       ES
       Vola
       DDRecovery
       DDSpeed
           
    end
    
    methods 
        
        function MDD_Distribution =  MDDDistribution(self, Process)
            
            %Compute the maximum drawdown distribution of the Process Input over a path
            % of length length(Path)

            %   Input: 

            % 1. Process : Timeseries of returns
            % 2. Path : Integer Number of days (all period are working)

            % Output: 

            % 1. Distribution : Distribution of Maximum Drawdown over the Path

            %******************************************************************

            % Number of days in the process
            [nDays,nAsset] = size(Process);

            % Transform Returns into price
            Price = ones(nDays+1,nAsset);
            Price(2:end,:) = cumprod(1+Process);

            % Preallocate Variable
            MDD_Distribution = zeros(nDays+1-self.PathLength,nAsset);

            % Computing the distribution
                for k = 1:nDays+1-self.PathLength
                MaxData = Price(k,:);
                MaxDD = zeros(1,nAsset);
                    for i = 1:self.PathLength
                    MaxData = max(MaxData, Price(k+i,:));
                    DD = (MaxData - Price(i+k,:)) ./ MaxData;
                        if any(DD > MaxDD)
                        p = DD > MaxDD;
                        MaxDD(p) = DD(p);
                        end
                    end
                MDD_Distribution(k,:) = MaxDD;
                end
        end
          
        
    end
    
    
end