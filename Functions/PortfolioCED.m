function [CED_Portfolio] = PortfolioCED(Weights, Data, Alpha, Window)

%Construction of Portfolio
Data = Weights*Data';
Dis = MDD_Distribution(Data', Window);
CED_Portfolio = zeros(size(Data,2)-2*Window, 1);
for i = 1:length(CED_Portfolio)
    [~, CED_Portfolio(i)] = CED(Dis(i:i+Window),Alpha);
end

end

