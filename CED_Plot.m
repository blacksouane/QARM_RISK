function [f,MDD,Quan,CondExpDD] = CED_Plot(Process,Alpha,Path)
%Plot of the distribution of DrawDOwn + CED
%   INPUT:

% 1. A return process
% 2. Alpha : Level of "Confidence" - Numeric Value between 0 ans 1
% 3. Path the Length of the path

% OUTPUT : 

%  f : a figure object


%Computing the distribution of DrawDown
MDD = MDD_Distribution(Process,Path);

%Computing the quantiles and CED
[Quan,CondExpDD] = CED(MDD,Alpha);

%Ploting the results 
f = figure('visible','on');
histogram(MDD,'Normalization','probability');
title(sprintf('Distribution of MDD with CED %s for path of length %s',string(Alpha*100) + '%',string(int16(Path))))
xlabel('MDD')
ylabel('Frequency')
xline(Quan,'-.b')
xline(CondExpDD,'-.r')
legend('Distribution of MDD',sprintf('%s quantile',string(Alpha)),sprintf('CED %s %',string(Alpha)))


end
