


clc 
close all


%% Data Import and preprocessing
ImportData;

Data.mkt = (table2array(FFResearchDataFactorsdaily(2:end,2))+ table2array(FFResearchDataFactorsdaily(2:end,5)))/100;
Data.rf = table2array(FFResearchDataFactorsdaily(2:end,5))/100;
Data.date = datetime(table2array(FFResearchDataFactorsdaily(2:end,1)),'ConvertFrom','yyyymmdd');

clear FFResearchDataFactorsdaily


%% Conditional Expected Drawdown

LenPath = 125;
Alpha = 0.05;
MDD = MDD_Distribution(Data.mkt,LenPath);
[Quan,CED_True] = CED(MDD,Alpha);


%SimulatedProcess Computation
Nspace = linspace(0,1,101);
Data.CED_SIMULATED = zeros(1,101);
Data.MDD_Simulated = zeros(25000-LenPath+1,101);
count = 1;

for i = Nspace
    
    SIMDATA = SimulatedProcess(25000,i);
    Data.MDD_Simulated(:,count) = MDD_Distribution(SIMDATA,LenPath);
    [~,Data.CED_SIMULATED(count)] = CED(Data.MDD_Simulated(:,count),Alpha);
    count = count + 1;
    
end 