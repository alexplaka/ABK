function [estimates] = FitCurve_Hill_fixK(xdata, ydata)
% Fits data to Hill equation y = x^n / (x^n + K)
% Parameter K is specified.
% Parameter n (the Hill coefficient, nH is optimized to fit the curve and 
% reported. A graph is produced.
% Input arguments xdata, ydata, must be row vectors.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

% Call fminsearch with a random starting point.
start_point = rand(1);
model = @fitfun;
estimates = fminsearch(model, start_point);
K = 31;
% fitfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for [xdata^n / (xdata^n + K)] - ydata, 
% and the FittedCurve. FMINSEARCH only needs sse, but we want to 
% plot the FittedCurve at the end.
    
    function [sse, FittedCurve] = fitfun(params)
        n = params;
        FittedCurve = 12 .* xdata.^n ./ (xdata.^n + 31);
        ErrorVector = FittedCurve - ydata;
        sse = sum(ErrorVector .^ 2);
    end

plot(xdata, ydata, '*')
hold on
[sse, FittedCurve] = model(estimates);
plot(xdata, FittedCurve, 'r')
 
xlabel('t')
ylabel('f(estimates,t)')
title(['Fitting PhosEvolution to Hill equation (K = ' num2str(K)]);
legend('data', ['y = 12 * x^n / (x^n + K)'],'Location','SouthEast');
xrange = xlim; xpos = xrange(2);
yrange = ylim; ypos = yrange(2)/4;

% Calculate Correlation Coefficient
Fit = 12 .* xdata.^estimates ./ (xdata.^estimates + K);
Compare = [(ydata)' (Fit)'];
[r p] = corrcoef(Compare);
R = r(1,2);     p_value = p(1,2);
% --------------------

text(xpos,ypos-1,['n=' num2str(estimates) ',K=' num2str(K)],...
    'HorizontalAlignment','right');
text(xpos,ypos,['R=' num2str(R)],'HorizontalAlignment','right');
hold off
end