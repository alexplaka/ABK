function [estimates,R] = FitCurve_lineb_s(xdata, ydata, y_int)
% Fits data to linear equation y = m * x + b (where b is input into the fc)
% Parameter m (the slope) is reported. Param. b (the y-intercept) is set.
% Input arguments xdata, ydata, MUST be *ROW* vectors.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

b = y_int;              % Set b to specified y-intercept

% Call fminsearch with a random starting point.
start_point = rand;
model = @fitfun;
estimates = fminsearch(model, start_point);

% fitfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for [m*xdata + b] - ydata, 
% and the FittedCurve. FMINSEARCH only needs sse, but we want to 
% plot the FittedCurve at the end.
    
    function [sse, FittedCurve] = fitfun(params)
        m = params(1);
        FittedCurve = m * xdata + b;
        ErrorVector = FittedCurve - ydata;
        sse = sum(ErrorVector .^ 2);
    end

[sse, FittedCurve] = model(estimates);
 
% Calculate Correlation Coefficient
Fit = estimates(1) * xdata + b;
Compare = [(ydata)' (Fit)'];
[r, p] = corrcoef(Compare);
R = r(1,2);     p_value = p(1,2)

end