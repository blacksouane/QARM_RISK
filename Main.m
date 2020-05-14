


clc 
close all
rng('default');

%import KevinShepperd Toolbox
addpath(genpath('C:\Users\Benjamin\OneDrive\1. HEC\Master\MScF 4.2\EMF\2020\Homeworks\KevinSheperdToolBox'));

% Import function library
addpath(genpath('C:\Users\Benjamin\OneDrive\Documents\GitHub\QARM_RISK'))

%% Data Import and preprocessing
disp('Import Data')

ImportData;

Data.mkt = (table2array(FFResearchDataFactorsdaily(2:end,2))+ table2array(FFResearchDataFactorsdaily(2:end,5)))/100;
Data.rf = table2array(FFResearchDataFactorsdaily(2:end,5))/100;
Data.date = datetime(table2array(FFResearchDataFactorsdaily(2:end,1)),'ConvertFrom','yyyymmdd');
Data.mkt = Data.mkt(1:end-2);
Data.rf = Data.rf(1:end-2);
Data.date = Data.date(1:end-2);
clear FFResearchDataFactorsdaily


%% Conditional Expected Drawdown

Alpha = 0.05;
LenPath = 125;

% ********************************************************************
% Path Length Relation
% ********************************************************************
disp('Path Length Relation')

Path.LenSpace = 100:10:1000;
count = 1;
Path.CED_Path = zeros(length(Path.LenSpace), 1);
for i = Path.LenSpace
    temp = MDD_Distribution(Data.mkt,i);
    [~, Path.CED_Path(count)] = CED(temp,Alpha);
    count = count + 1;
end
clear count

Path.PathRelation = fitlm(Path.LenSpace', Path.CED_Path, 'VarNames', {'Path Length', 'CED'});

f = figure('visible', 'on');
plot(Path.LenSpace, Path.CED_Path)
xlabel('Path Length')
ylabel('CED')
title('CED for different Path Length')
legend('CED at 5%')
print(f, 'Plots/PathLength','-dpng','-r1000')


% TAKE A LOT OF TIME, ONLY USEFUL TO DO A 3D PLOT
% Path.LenSpace = 100:10:5000;
%  Path.AlphaSpace = 0.01:0.01:0.99;
%  Path.PathAlphaMat = zeros(length(Path.LenSpace), length(Path.AlphaSpace));
%  countI = 1;
%  countJ = 1;
%  for i = Path.LenSpace
%     disp(i)
%     temp = MDD_Distribution(Data.mkt,i);
%     for j = Path.AlphaSpace
%     [~, Path.PathAlphaMat(countI, countJ)] = CED(temp,j);
%     countJ = countJ + 1;
%     end
%      countI = countI + 1;
%      countJ = 1;
% end
% 
% f = figure('visible','on');
% colormap(hot(256));
% surf(Path.AlphaSpace, Path.LenSpace, Path.PathAlphaMat)
% camlight right;
% lighting phong;
% shading interp
% xlabel('1 - Alpha')
% ylabel('Number of Day')
% zlabel('CED')
% print(f, 'Plots/3D_AlphaPathCED','-dpng', '-r1000')

% ********************************************************************
% Frequency of observations
% ********************************************************************
disp('Frequency Relation')

Frequency.FreqSpace = 100:10:1000; 
Frequency.CED_Freq = zeros(length(Frequency.FreqSpace),3);
count = 1;
for i = Frequency.FreqSpace
    
    for j = 1:3
        if j == 1
            temp = MDD_Distribution(Data.mkt,i);
            [~, Frequency.CED_Freq(count, j)] = CED(temp,Alpha); 
        elseif j == 2
            tempData = FrequencyConverter(Data.mkt, j-1);
            temp = MDD_Distribution(tempData,round(i/5));
            [~, Frequency.CED_Freq(count, j)] = CED(temp,Alpha); 
        else
            tempData = FrequencyConverter(Data.mkt, j-1);
            temp = MDD_Distribution(tempData,round(i/20));
            [~, Frequency.CED_Freq(count, j)] = CED(temp,Alpha); 
        end 
    end 
    count = count + 1;
end

%Reshape the data and perform a linear regression
Frequency.toFit = reshape(Frequency.CED_Freq,[length(Frequency.FreqSpace)*3, 1]);
Frequency.Freq = zeros([length(Frequency.FreqSpace)*3, 2]);
for i = 1:3
   Frequency.Freq((i-1)*length(Frequency.FreqSpace)+1:i*length(Frequency.FreqSpace), 1)...
       = Frequency.FreqSpace';
   Frequency.Freq((i-1)*length(Frequency.FreqSpace)+1:i*length(Frequency.FreqSpace), 2)...
       = i;
end

Frequency.FrequencyRelation = fitlm(Frequency.Freq, Frequency.toFit, 'CategoricalVars', logical([0, 1]), ...
    'VarNames', {'Path Length', 'Frequency', 'CED'});

f = figure('visible', 'on');
plot(Frequency.FreqSpace, Frequency.CED_Freq(:, 1))
hold on
plot(Frequency.FreqSpace, Frequency.CED_Freq(:, 2))
hold on
plot(Frequency.FreqSpace, Frequency.CED_Freq(:, 3))
xlabel('Path Length')
ylabel('CED')
title('CED at 5% for different Path Length')
legend('Daily Frequency', 'Weekly Frequency', 'Monthly Frequency','location','best')
print(f, 'Plots/Frequency','-dpng','-r1000')

% ********************************************************************
% Auto-Correlations
% ********************************************************************
disp('Simulated AR1 Process')

%SimulatedProcess Computation
Nspace = linspace(0,1,101);
Data.CED_SIMULATED = zeros(1,101);
Data.Vola_SIMULATED = zeros(1,101);
Data.ES_SIMULATED = zeros(1,101);
Data.MDD_Simulated = zeros(15000-LenPath+1,101);
count = 1;

for i = Nspace
    
    SIMDATA = SimulatedProcess(15000,i);
    Data.MDD_Simulated(:,count) = MDD_Distribution(SIMDATA,LenPath);
    [~,Data.CED_SIMULATED(count)] = CED(Data.MDD_Simulated(:,count),Alpha);
    Data.Vola_SIMULATED(count) = std(SIMDATA);
    [~, Data.ES_SIMULATED(count), ~] = ES(SIMDATA,Alpha);
    count = count + 1;
    
end 
clear count

% Regression
ArRelation = fitlm(Nspace', Data.CED_SIMULATED','VarNames', {'Kappa', 'CED'});

f = figure('visible', 'on');
yyaxis left
ylabel('CED Scale')
plot(Nspace(1:end-30), smooth(Data.CED_SIMULATED(1:end-30)))
hold on
yyaxis right
plot(Nspace(1:end-30), smooth(Data.ES_SIMULATED(1:end-30)))
hold on
plot(Nspace(1:end-30), smooth(Data.Vola_SIMULATED(1:end-30)))
xlabel('AutoCorrelation Factor')
ylabel('Vol and ES scale')
title('Measure of risk for different level of AutoCorrelation')
legend('CED', 'ES', 'Vol','location','best')
print(f, 'Plots/AutoCorrelation','-dpng','-r1000')

% ********************************************************************
% AR1 fit rolling window
% ********************************************************************
AR.LenPath = 125;
AR.Alpha = 0.05;
AR.NumFit = length(Data.mkt) - AR.LenPath;
AR.ARfit = zeros(AR.NumFit, 1);

for i = 1:AR.NumFit
    temp = armaxfilter(Data.mkt(i:i+AR.LenPath),1,1);
    AR.ARfit(i) = temp(2);
end

AR.Dis = MDD_Distribution(Data.mkt, AR.LenPath);
AR.CED = zeros(AR.NumFit-AR.LenPath, 1);
for i = 1:length(AR.ARfit)-AR.LenPath
    [~, AR.CED(i)] = CED(AR.Dis(i:i+AR.LenPath),AR.Alpha);
end

% k - dd Space
f = figure('visible','on');
scatter(AR.CED(20000:end),AR.ARfit(20000+AR.LenPath:end))
title('AutoCorrelation and CED relation')
ylabel('Kappa (k)')
xlabel('CED')
legend('k - CED space','location','best')
print(f, 'Plots/AR_CED_RolWindow','-dpng','-r1000')

f = figure('visible','on');
plot(AR.Dis(1:end-1))
hold on
plot(AR.ARfit)
title('AutoCorrelation and CED relation')
ylabel('Kappa (k)')
xlabel('CED')
legend('CED space','Kappa','location','best')
print(f, 'Plots/AR_CED_RolWindowDistribution','-dpng','-r1000')


%% Peak,speed and Recovery 
disp('Peak, Speed and Recovery Relation')

[MaxDD.MDD, MaxDD.Idx, MaxDD.Recover, MaxDD.SpeedRecover, MaxDD.SpeedMDD] = MaxDD_Speed(Data.mkt,1000);

% Distribution of the speeds
f = figure('visible','on');
histogram(MaxDD.SpeedRecover,100,'Normalization','probability');
hold on
histogram(MaxDD.SpeedMDD,25,'Normalization','probability');
title('Speed of MDD and Recovery')
xlabel('Ndays')
ylabel('Frequency')
legend('Recovery Speed', 'Drawdown Speed','location','best')
print(f, 'Plots/SpeedRecoveryHist','-dpng','-r1000')

% Speed vector
f = figure('visible','on');
plot(Data.date(1001:end),MaxDD.SpeedRecover);
hold on
plot(Data.date(1001:end),MaxDD.SpeedMDD);
title('Speed of MDD and Recovery')
xlabel('Date')
ylabel('Ndays')
legend('Recovery Speed', 'Drawdown Speed','location','best')
print(f, 'Plots/SpeedRecoveryLine','-dpng','-r1000')

% Analysing the relation
Prediction.LM = fitlm([MaxDD.MDD, MaxDD.SpeedMDD], MaxDD.SpeedRecover,...
    'VarNames', {'Intensity', 'Speed','Recovery'});

% Gaussian Process
Prediction.Split = 0.8;
Prediction.X = [MaxDD.MDD, MaxDD.SpeedMDD];
Prediction.Mean = mean(Prediction.X);
Prediction.Std = std(Prediction.X);
Prediction.X = (Prediction.X - Prediction.Mean)./Prediction.Std;
Prediction.Indices = randperm(size(MaxDD.Recover,1));
Prediction.MeanY = mean(rmmissing(MaxDD.SpeedRecover));
Prediction.StdY = std(rmmissing(MaxDD.SpeedRecover));
Prediction.y = (MaxDD.SpeedRecover-Prediction.MeanY)./Prediction.StdY;
Prediction.yTrain = Prediction.y(Prediction.Indices(1:round(Prediction.Split...
    *size(MaxDD.Recover,1))));
Prediction.yTest = Prediction.y(Prediction.Indices(round(Prediction.Split...
    *size(MaxDD.Recover,1))+1:end));
Prediction.XTrain = Prediction.X(Prediction.Indices(1:round(Prediction.Split...
    *size(MaxDD.Recover,1))))';
Prediction.XTest = Prediction.X(Prediction.Indices(round(Prediction.Split...
    *size(MaxDD.Recover,1))+1:end))';
Prediction.GPR = fitrgp(Prediction.XTrain,Prediction.yTrain,...
    'BasisFunction','pureQuadratic',...
    'Standardize',true,'KernelFunction','matern32');
Prediction.yPred = predict(Prediction.GPR,Prediction.XTest);
Prediction.Loss = loss(Prediction.GPR, Prediction.XTest,Prediction.yTest);

% Plot of the prediction
f = figure('visible','on');
scatter(Prediction.yPred,Prediction.yTest)
hold on
plot([0 5000], [0 5000])
xlabel('Prediction')
ylabel('Actual Value')
xlim([0 5000]) % Reduce to scale to exclude high Drawdown 
ylim([0 5000]) % Reduce to scale to exclude high Drawdown 
title('Recovery Period prediction')
legend('Prediction X True value','Perfect prediction Line','location','southeast')
print(f, 'Plots/SpeedRecoveryPred','-dpng','-r1000')

% Plot of the prediction
f = figure('visible','on');
plot(Prediction.yTest,'.')
hold on
plot(Prediction.yPred,'.')
xlabel('Data Points')
ylabel('Recovery Speed')
title('Recovery Period prediction')
legend('True Value','Prediction','location','southeast')
print(f, 'Plots/SpeedRecoveryPred','-dpng','-r1000')

clear i j count SIMDATA temp tempData f

% 3d plot
xNodes = 0:0.001:1;
yNodes = linspace(10, 850, 1001);
z = gridfit(MaxDD.MDD, MaxDD.SpeedMDD, MaxDD.SpeedRecover, xNodes, yNodes,...
    'smoothness', 0.7);

figure = figure('visible','on');
colormap(hot(256));
surf(xNodes,yNodes,z);
camlight right;
lighting phong;
shading interp
xlabel('DrawDown in %')
ylabel('Speed of the DrawDown in Days')
zlabel('Recovery Period in Days')
title 'Recovery Speed w.r.t to the speed and intensity of the DD'
print(figure, 'Plots/3d','-dpng','-r1000')
%% Risk Contribution

ImportIndustryData;
DateImport;
Industry.Data = table2array(IndustryPortfoliosDaily(1:end-1,2:end))./100;
toConvert = table2array(IndustryPortfoliosDaily_2(2:end,1));
toConvert = erase(toConvert,',');
Industry.Date = datetime(toConvert,'InputFormat','yyyyMMdd');
clear IndustryPortfoliosDaily IndustryPortfoliosDaily_2 toConvert
Industry.Data = Industry.Data(1:24704,:);
Industry.Date = Industry.Date(1:24704);
Industry.Names = {'Consumer NonDurables', 'Consumer Durables', 'Manufacturing', ...
    'Energy, Oil, Gas', 'Tech', 'Telecom', 'Wholesale and Retail', ...
    'Healthcare', 'Utilities','Other'};
% Equally weighted Portfolio
[Industry.NumDays, Industry.NumAsset] = size(Industry.Data);
Industry.Alpha = 0.05;
Industry.LenPath = 125;
Industry.EQ = sum(Industry.Data,2)./Industry.NumAsset;

% Marginal Risk Contribution
h = 0.001;
Window = 500;
Industry.EqWeights = ones(1,Industry.NumAsset)./Industry.NumAsset;
[Industry.MCR, Industry.CED] = RiskContribution(Industry.Data,...
    Industry.EqWeights, Industry.Alpha, Window,h);

Industry.FRC = Industry.MCR./10;
for days = 1:length(Industry.MCR)
    Industry.FRC(days, :) = Industry.FRC(days, :)./Industry.CED(days);
    Industry.FRC(days, :) = Industry.FRC(days, :)./sum(Industry.FRC(days, :));
end


f = figure('visible', 'on');
area(Industry.Date(1001:end),Industry.FRC)
ylim([-0.1 1.1])
legend(Industry.Names, 'location','bestoutside')
title('Fractional Contribution to risk of each industry')
print(f, 'Plots/MCRarea','-dpng','-r1000')
