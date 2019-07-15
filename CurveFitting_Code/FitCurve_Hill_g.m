function [estimates,Fit] = FitCurve_Hill_g(xdata, ydata, maxy)
% Fits data to Hill equation y = A * x^n / (x^n + K)
% Parameters n (the Hill coefficient, nH) and K are reported.
% Input arguments xdata, ydata, must be ROW vectors.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

if nargin == 3
    A = maxy;
elseif nargin == 2
    A = max(ydata);           % Scale fitted curve to data maximum
end

% Call fminsearch with a random starting point.
start_point = rand(2);
model = @fitfun;
estimates = fminsearch(model, start_point);

% fitfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for [xdata^n / (xdata^n + K)] - ydata, 
% and the FittedCurve. FMINSEARCH only needs sse, but we want to 
% plot the FittedCurve at the end.
    
    function [sse, FittedCurve] = fitfun(params)
        n = params(1);
        K = params(2);
        FittedCurve = A * xdata.^n ./ (xdata.^n + K);
        ErrorVector = FittedCurve - ydata;
        sse = sum(ErrorVector .^ 2);
    end

% ---------- Draw Graph -----------
% figure;
semilogx(xdata, ydata, '.','MarkerSize',5.0)
% plot(xdata, ydata, '.','MarkerSize',5.0)
axis([min(xdata) max(xdata) 0 1]);
hold on;
[sse, FittedCurve] = model(estimates);
plot(xdata, FittedCurve, 'r')
 
% xlabel('t');
% ylabel('f(estimates,t)');
ylabel('Fraction of subunits');
title('Fit of ---- to Hill equation');
legend('data', 'y = A * x^n / (x^n + K)','Location','SouthEast');
xrange = xlim;                      xpos = xrange(2);
yrange = ylim;                      ypos = yrange(2)/4;
% ---------------------------------

% Calculate Correlation Coefficient
Fit = A * xdata.^estimates(1) ./ (xdata.^estimates(1) + estimates(2));
Compare = [(ydata)' (Fit)'];
[r p] = corrcoef(Compare);
R = r(1,2);     p_value = p(1,2);
% Calculate sensitivity coefficient (Rj)
% Ratio of concentration of effector needed to produce a 90% response
% to that producing a 10% response.
x90 = (9 * estimates(2))^(1/estimates(1));
x10 = (1/9 * estimates(2))^(1/estimates(1));
Rj = x90 / x10;

% --------- Finish graph -----------
text(xpos,ypos,{['A=' num2str(A)] ...
    ['nH=' num2str(estimates(1)) ',K=' num2str(estimates(2))] ...
    ['R_j=' num2str(Rj)]}, 'HorizontalAlignment','right');
text(xpos,ypos+0.1,['R=' num2str(R)],'HorizontalAlignment','right');
hold off;

% Print results to std output
disp(' ');
disp(['A=' num2str(A)]);
disp(['nH=' num2str(estimates(1)) ', K=' num2str(estimates(2)) ...
    ', R_j=' num2str(Rj)]);
disp(['R=' num2str(R) ' p=' num2str(p_value)]); 
disp(' ');

end