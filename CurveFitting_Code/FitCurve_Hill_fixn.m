function [estimates,FittedCurve] = FitCurve_Hill_fixn(xdata, ydata)
% Fits data to Hill equation y = x^n / (x^n + K)
% Parameters n (the Hill coefficient, nH) is specified.
% Parameter K is optimized to fit the curve and reported.
% Input arguments xdata, ydata, must be row vectors.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

n = 1;      % Hill coefficient, nH, is set to this value

% Call fminsearch with a random starting point.
start_point = rand(1);
model = @fitfun;
estimates = fminsearch(model, start_point);

% fitfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for [xdata^n / (xdata^n + K)] - ydata, 
% and the FittedCurve. FMINSEARCH only needs sse, but we want to 
% plot the FittedCurve at the end.
    
    function [sse, FittedCurve] = fitfun(params)
        K = params;
        FittedCurve = xdata.^n ./ (xdata.^n + K);
        ErrorVector = FittedCurve - ydata;
        sse = sum(ErrorVector .^ 2);
    end

% ---------- Plot graph ------------
figure;
% semilogx(xdata, ydata, '.','MarkerSize',4.0)
plot(xdata, ydata, '.','MarkerSize',4.0);
axis tight;
hold on;
[sse, FittedCurve] = model(estimates);
plot(xdata, FittedCurve, 'r');
 
% xlabel('t')
% ylabel('f(estimates,t)')
title(['Fit of ---- to Hill equation (n = ' num2str(n) ')']);
legend('data', ['y = x^' num2str(n) ' / (x^' num2str(n) ' + K)'],'Location','SouthEast');
xrange = xlim; xpos = xrange(2);
yrange = ylim; ypos = yrange(2)/4;

% Calculate Correlation Coefficient
Fit = xdata.^n ./ (xdata.^n + estimates);
Compare = [(ydata)' (Fit)'];
[r p] = corrcoef(Compare);
R = r(1,2);     p_value = p(1,2);
% -------------------------------------------------
% Calculate sensitivity coefficient (Rj)
% Ratio of concentration of effector needed to produce a 90% response
% to that producing a 10% response.
x90 = (9 * estimates)^(1/n);
x10 = (1/9 * estimates)^(1/n);
Rj = x90 / x10;

text(xpos,ypos,['K=' num2str(estimates)],'HorizontalAlignment','right');
text(xpos,ypos+0.1,['R=' num2str(R)],'HorizontalAlignment','right');
hold off;

% Print results to std output
disp(' ');
disp(['nH=' num2str(n) ', K=' num2str(estimates) ...
    ', R_j=' num2str(Rj)]);
disp(['R=' num2str(R) ', p=' num2str(p_value)]); 

end