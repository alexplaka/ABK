function [estimates] = FitCurve_sigmoid(xdata, ydata)
% Fits data to equation y = A (1 - e^(-kx))^h
% Parameters k and h are optimized to fit the curve and reported.
% Input arguments xdata, ydata, must be COLUMN vectors.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

A = max(ydata)           % Scale fitted curve to data maximum
% Call fminsearch with a random starting point.
start_point = rand(1, 2);
model = @fitfun;
estimates = fminsearch(model, start_point);
% fitfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for A * exp(-lambda * xdata) - ydata, 
% and the FittedCurve. FMINSEARCH only needs sse, but we want to 
% plot the FittedCurve at the end.
    function [sse, FittedCurve] = fitfun(params)
        k = params(1);
        h = params(2);
        FittedCurve = A .* (1 - exp(-k .* xdata)).^h;
        ErrorVector = FittedCurve - ydata;
        sse = sum(ErrorVector .^ 2);
    end

[sse, FittedCurve] = model(estimates);
% Calculate Correlation Coefficient
Fit = A .* (1 - exp(-estimates(1) .* xdata)).^h;
Compare = [ydata Fit];
[r p] = corrcoef(Compare);
R = r(1,2);     p_value = p(1,2);

% Plot
semilogx(xdata, ydata, '*');                    hold on;
semilogx(xdata, FittedCurve, 'r'); 
xlabel('x');                                ylabel('f(estimates,x)')
% title(['Fitting to function (1 - exp^(-k * x))^h']);
legend('data', ['y = A (1 - e^-^k^x)^h'],'Location','Best');
xrange = xlim; xpos = xrange(2);
yrange = ylim; ypos = 3*yrange(2)/4;
text(xpos,ypos,['k=' num2str(estimates(1)) ' ; h=' num2str(estimates(2))...
    ' R=' num2str(R) ', p=' num2str(p_value)], ...
    'HorizontalAlignment','right');
hold off

end