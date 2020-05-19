function [CED_Portfolio] = PortfolioCEDGeneral(Weights, Data, Alpha, Window)

%Construction of Portfolio
Data = Weights*Data';
CED_Portfolio = CED_Faster(Data', Window, Alpha);
end

