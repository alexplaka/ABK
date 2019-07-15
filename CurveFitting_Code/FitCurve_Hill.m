function [estimates,Rj] = FitCurve_Hill(xdata, ydata)
% This function fits the data to the Hill equation and the results are
% printed on standard output.  No graphs of the data are plotted.
% Parameters n (the Hill coefficient, nH) and K are reported.
% A graph is NOT produced.
% Input arguments xdata, ydata, must be row vectors.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

% Call fminsearch with a random starting point.
start_point = rand(2);
model = @fitfun;
estimates = fminsearch(model, start_point);

% fitfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for [xdata^n / (xdata^n + K)] - ydata, 
% and the FittedCurve. FMINSEARCH only needs sse.
    
    function [sse, FittedCurve] = fitfun(params)
        n = params(1);
        K = params(2);
        FittedCurve = xdata.^n ./ (xdata.^n + K);
        ErrorVector = FittedCurve - ydata;
        sse = sum(ErrorVector .^ 2);
    end

% Calculate Correlation Coefficient
Fit = xdata.^estimates(1) ./ (xdata.^estimates(1) + estimates(2));
Compare = [(ydata)' (Fit)'];
[r p] = corrcoef(Compare);
R = r(1,2);     p_value = p(1,2);
% Calculate sensitivity coefficient (Rj)
% Ratio of concentration of effector needed to produce a 90% response
% to that producing a 10% response.
x90 = (9 * estimates(2))^(1/estimates(1));
x10 = (1/9 * estimates(2))^(1/estimates(1));
Rj = x90 / x10;

% Print results to std output
disp(' ');
disp(['nH=' num2str(estimates(1)) ', K=' num2str(estimates(2)) ...
    ', R_j=' num2str(Rj)]);
disp(['R=' num2str(R) ', p=' num2str(p_value)]); 

end