function [estimates] = FitCurve_line_s(xdata, ydata)
% Fits data to linear equation y = m * x + b
% Parameters m (the slope) and b (the y-intercept) are reported.
% Input arguments xdata, ydata, MUST be *ROW* vectors.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

% Call fminsearch with a random starting point.
start_point = rand(2);
model = @fitfun;
estimates = fminsearch(model, start_point);

% fitfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for [xdata^n / (xdata^n + K)] - ydata, 
% and the FittedCurve. FMINSEARCH only needs sse, but we want to 
% plot the FittedCurve at the end.
    
    function [sse, FittedCurve] = fitfun(params)
        m = params(1);
        b = params(2);
        FittedCurve = m * xdata + b;
        ErrorVector = FittedCurve - ydata;
        sse = sum(ErrorVector .^ 2);
    end

% Calculate Correlation Coefficient
Fit = estimates(1) * xdata + estimates(2);
Compare = [(ydata)' (Fit)'];
[r p] = corrcoef(Compare);
R = r(1,2);     p_value = p(1,2);

% ---------- Draw Graph -----------
% figure;
% plot(xdaRta, ydata, '.','MarkerSize',5.0);                   hold on;

% [sse, FittedCurve] = model(estimates);
% plot(xdata, FittedCurve, 'r');

% axis([min(xdata) max(xdata) 0 1]);                    
% xlabel('t');
% ylabel('Fraction of subunits');
% legend('data', 'y = m * x + b','Location','SouthEast');
% xrange = xlim;                      xpos = xrange(2);
% yrange = ylim;                      ypos = yrange(2)/4;

% text(xpos-0.1,ypos,['m=' num2str(estimates(1)) ',b=' num2str(estimates(2))]);
% text(xpos,ypos+0.1,['R=' num2str(R)],'HorizontalAlignment','right');
% hold off;
% ---------------------------------

% Print results to std output
disp(' ');
disp(['m=' num2str(estimates(1)) ', b=' num2str(estimates(2))]);
disp(['R=' num2str(R) ' p=' num2str(p_value)]); 
disp(' ');

end