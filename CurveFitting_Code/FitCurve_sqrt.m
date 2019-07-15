function [estimates] = FitCurve_sqrt(xdata, ydata)
% Fits data to equation y = k x^(1/2)
% Parameter k is optized to fit the data and reported.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

% Call fminsearch with a random starting point.
start_point = rand(1);
model = @fitfun;
estimates = fminsearch(model, start_point);
% fitfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for [k * xdata^0.5] - ydata, 
% and the FittedCurve. FMINSEARCH only needs sse, but we want to 
% plot the FittedCurve at the end.
    
    function [sse, FittedCurve] = fitfun(params)
        k = params;
        FittedCurve = k * sqrt(xdata);
        ErrorVector = FittedCurve - ydata;
        sse = sum(ErrorVector .^ 2);
    end

% plot(xdata, ydata, '*')
hold on
[sse, FittedCurve] = model(estimates);
plot(xdata, FittedCurve, 'r')
 
xlabel('x')
ylabel('f(estimates,x)')
% title(['Fitting to function k * sqrt(x)']);
legend('data', ['y = k * sqrt(x)']);
xrange = xlim; xpos = xrange(2);
yrange = ylim; ypos = 3*yrange(2)/4;

% Calculate Correlation Coefficient
Fit = estimates * sqrt(xdata);
Compare = [ydata Fit];
[r p] = corrcoef(Compare)
R = r(1,2);
% --------------------

text(xpos,ypos,['k=' num2str(estimates) ', R=' num2str(R)],'HorizontalAlignment','right');
hold off
end