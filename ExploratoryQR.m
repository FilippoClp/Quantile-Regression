%% Clear workspace, read data from .xls file and insert quantile and quarter parameters

clear all
close all
clc

% Read excel file and extract variables
[Data,VarNames] = xlsread('data.xlsx', 'Data');

%Time = datenum(TEXT(2:end, 1),'dd/mm/yyyy');
VarNames = VarNames(1,2:end); 

%Set number of quarters ahead
h=4;

%Forecast settings 
QLow=0.1;
QMid=0.5;
QHigh=0.9;

%Number of Bootstrap replications
Boot=500;
   
%% Quantile regression estimates

%Initialize table and create Row names
CoeffName = {'OLS';strcat('Q',num2str(QLow*100));strcat('Q',num2str(QMid*100));strcat('Q',num2str(QHigh*100)); strcat('Tick_loss_Q',num2str(QLow*100))}; 
RowName = {}; OLSest = {}; Q10est = {}; Q50est = {}; Q90est = {}; TickLoss = {};

%Loop through the regressors and build the table
for i = 2:length(VarNames) 
    
    %Create dependant variable and regressors
    Y = Data((h + 1):end , 1); 
    X = Data(1:(end - h), i);
    
    %OLS regression
    modelOLS = fitlm(X, Y); 
      
    %Quantile regressions
    [p10,stats10] = quantreg(X, Y, QLow, 1, Boot);
    [p50,stats50] = quantreg(X, Y, QMid, 1, Boot);
    [p90,stats90] = quantreg(X, Y, QHigh, 1, Boot);
    
    %Alternative calculation with rq, which yields the same results!
    %Z=[ones(size(X)), X];
    %regr = rq(Z, Y, 0.10);
    
    %Calculate Tick Loss
    fit_q = polyval(p10, X);
    fcast_err = Y - fit_q;
    rho=@(r)abs(r.*(QLow -(r<0)));
    TL = mean(rho(fcast_err));
    
    %Calculate T tests and p-values for the quantile regression coefficient
    t = (p10(1))/(stats10.pse(1)); pvalue10 = 2*tcdf(t,length(X)); 
    t = (p50(1))/(stats50.pse(1)); pvalue50 = 2*tcdf(t,length(X));
    t = (p90(1))/(stats90.pse(1)); pvalue90 = 2*tcdf(t,length(X)); 
      
    %Table Coefficients of OLS
    OLSco = num2str(modelOLS.Coefficients.Estimate(2),5);
        if (abs(modelOLS.Coefficients.pValue(2)) < 0.1) 
            OLSco = strcat(OLSco,'*');
            if (abs(modelOLS.Coefficients.pValue(2)) < 0.05) 
                OLSco = strcat(OLSco,'*');
                if (abs(modelOLS.Coefficients.pValue(2)) < 0.001) 
                    OLSco = strcat(OLSco,'*'); 
                end 
            end 
        end
    
    %Table Coefficients of p10
    Q10co = num2str(p10(1),5);
        if (abs(pvalue10) < 0.1) 
            Q10co = strcat(Q10co,'*'); 
            if (abs(pvalue10) < 0.05) 
                Q10co = strcat(Q10co,'*'); 
                if (abs(pvalue10) < 0.01) 
                    Q10co = strcat(Q10co,'*'); 
                end 
            end 
        end
    
    %Table Coefficients of p50
    Q50co = num2str(p50(1),5);
        if (abs(pvalue50) < 0.1) 
            Q50co = strcat(Q50co,'*'); 
            if (abs(pvalue50) < 0.05) 
                Q50co = strcat(Q50co,'*'); 
                if (abs(pvalue50) < 0.01) 
                    Q50co = strcat(Q50co,'*'); 
                end 
            end 
        end
    
    %Table Coefficients of p90
    Q90co = num2str(p90(1),5);
        if (abs(pvalue90) < 0.1) 
            Q90co = strcat(Q90co,'*'); 
            if (abs(pvalue90) < 0.05) 
                Q90co = strcat(Q90co,'*'); 
                if (abs(pvalue90) < 0.01) 
                    Q90co = strcat(Q90co,'*'); 
                end 
            end 
        end

    %Standard errors of OLS, p10, p50 & p90
    OLSse = strcat('(',num2str(abs(modelOLS.Coefficients.SE(2)),5),')');
    Q10se = strcat('(',num2str(abs(stats10.pse(1)),5),')');
    Q50se = strcat('(',num2str(abs(stats50.pse(1)),5),')');
    Q90se = strcat('(',num2str(abs(stats90.pse(1)),5),')');

    %Concatenate arrays for table
    RowName = [RowName; VarNames(i); strcat('SE-',num2str(i))];
    OLSest = [OLSest;{OLSco; OLSse}];
    Q10est = [Q10est;{Q10co; Q10se}];
    Q50est = [Q50est;{Q50co; Q50se}];
    Q90est = [Q90est;{Q90co; Q90se}];
    TickLoss = [TickLoss;{num2str(TL); num2str(NaN)}];
    
    %Plot Scatterplot
    f = figure('visible','off');
    scatter(X,Y); box on; hold on;
    plot(X, modelOLS.Coefficients.Estimate(1) + modelOLS.Coefficients.Estimate(2)*X, ...
        X, p10(2) + p10(1)*X, ...
        X, p50(2) + p50(1)*X, ...
        X, p90(2) + p90(1)*X);
    ylabel(VarNames(1)); xlabel(VarNames(i)); 
    legend('Scatter','OLS', strcat('Q',num2str(QLow*100)), strcat('Q',num2str(QMid*100)), strcat('Q',num2str(QHigh*100))); 
    axis tight; hold off;
    saveas(f,fullfile('Scatterplots',string(strcat(VarNames(i),'.png'))));
    
    %Plot Vulnerability Band
    %Res = QRboot(Data(:, i), Data(:, 1), [1, 4], 1, [0.10:0.10:0.90], Boot, 4);
    %BQ  = Res.BQ(2, :, h);
    %B2  = Res.B2(2, h);
    %bBQ = squeeze(Res.bBQ(2, :, h, :));
    %bB2 = squeeze(Res.bB2(2, h, :));
    %QQ = Res.QQ;
    %filename = fullfile('Nonlinearity Tests', string(strcat(VarNames(i), '_Q', num2str(h), '.png')));
    %PlotQRbands(BQ, bBQ, B2, bB2, QQ, 'filename', filename);   
         
end

TableFull = table(OLSest,Q10est,Q50est,Q90est,TickLoss,'RowNames',RowName,'VariableNames',CoeffName)
writetable(TableFull,'Table.xlsx','WriteRowNames',true);