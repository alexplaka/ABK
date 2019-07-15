function [estimates] = FitCurve_lineb_g(xdata, ydata, y_int)
% Fits data to linear equation y = m * x + b (where b is input into the fc)
% Parameters m (the slope) and b (the y-intercept) are reported.
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

% ---------- Draw Graph -----------
% figure;
plot(xdata, ydata, '.','MarkerSize',5.0);                   hold on;
% axis([min(xdata) max(xdata) 0 1]);                    
[sse, FittedCurve] = model(estimates);
plot(xdata, FittedCurve, 'r');
 
% xlabel('t');
ylabel('Fraction of subunits');
legend('data', 'y = m * x + b','Location','SouthEast');
xrange = xlim;                      xpos = xrange(2);
yrange = ylim;                      ypos = yrange(2)/4;
% ---------------------------------

% Calculate Correlation Coefficient
Fit = estimates(1) * xdata + b;
Compare = [(ydata)' (Fit)'];
[r p] = corrcoef(Compare);
R = r(1,2);     p_value = p(1,2);

% --------- Finish graph -----------
text(xpos-0.1,ypos,['m=' num2str(estimates(1)) ',b=' num2str(b)]);
text(xpos,ypos+0.1,['R=' num2str(R)],'HorizontalAlignment','right');
hold off;

% Print results to std output
disp(' ');
disp(['m=' num2str(estimates(1)) ', b=' num2str(b)]);
disp(['R=' num2str(R) ' p=' num2str(p_value)]); 
disp(' ');

end