function [estimates] = FitCurve_rec(xdata, ydata)
% Fits data to equation y = 1 / (kx)
% Parameter k is optimized to fit the curve and reported.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

% Call fminsearch with a random starting point.
start_point = rand(1);
model = @fitfun;
estimates = fminsearch(model, start_point);
% fitfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for [1 / (k * xdata)] - ydata, 
% and the FittedCurve. FMINSEARCH only needs sse, but we want to 
% plot the FittedCurve at the end.
    function [sse, FittedCurve] = fitfun(params)
        k = params;
        FittedCurve = 1 ./ (k * xdata);
        ErrorVector = FittedCurve - ydata;
        sse = sum(ErrorVector .^ 2);
    end

plot(xdata, ydata, '*')
hold on
[sse, FittedCurve] = model(estimates);
plot(xdata, FittedCurve, 'r')
 
xlabel('x')
ylabel('f(estimates,x)')
% title(['Fitting to function 1 / (k * x)']);
legend('data', ['y = 1 / (k * x)']);
xrange = xlim; xpos = xrange(2);
yrange = ylim; ypos = 3*yrange(2)/4;

% Calculate Correlation Coefficient
Fit = 1 ./ (estimates * xdata);
Compare = [ydata Fit];
[r p] = corrcoef(Compare)
R = r(1,2);
% --------------------

text(xpos,ypos,['k=' num2str(estimates) ', R=' num2str(R)],'HorizontalAlignment','right');
hold off
end