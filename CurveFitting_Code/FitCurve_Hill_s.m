function [estimates,Fit] = FitCurve_Hill_s(xdata, ydata)
% Fits data to Hill equation y = x^n / (x^n + K)
% Parameters n (the Hill coefficient, nH) and K are reported.
% This function differs from 'FitCurve_Hill.m' in that it has no output to
% the command console. 's' stands for silent. Also, it produces no graph of
% the fitted data.
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

end