% LRFIT_param_conf
% Alexis code, 16/08/2016
% coeffs(1) is slope, coeffs(2) is saturation and coeffs(3) is Rd.
% confint_95(1,1) is lower interval of slope

function [coeffs,confint_95] = LRFIT_param_confint(x,y,condition)

%% Fit: 'Winter'.    
[xData, yData] = prepareCurveData( x(condition), y(condition) );
if length(xData(xData<500)) < 8 % enough data close to PAR = 0 
    coeffs = NaN;
    confint_95 = NaN;
else
    % Set up fittype and options.
    % Mitscherlich
    ft = fittype( '-(b+c)*(1-exp((-a*x)/(b+c)))+c', 'independent', 'x', 'dependent', 'y' );
    % Michaelis-Menten
    %ft = fittype( '(a*b*x)/(a*x+b)+c', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [-Inf -Inf 0.2];
    opts.Upper = [Inf Inf Inf];
    opts.StartPoint = [0.01 6 2];

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    coeffs = coeffvalues(fitresult);
    confint_95 = confint(fitresult);
end

% Plot fit with data.
%figure( 'Name', 'Winter' );
% ax = gca;
% ax.FontSize = 16;
% ax.XLim = [0 3000];
% ax.YLim = [-20 15];
% h = plot( fitresult,'k');% xData, yData );

end