function [estimates, FittedCurve] = FitCurve_exp_Fixk(xdata, ydata, ko)
% Fits data to equation y = A e^(-kx)
% Produces NO graphical output.
% Parameter A are optimized to fit the curve and reported.
% Input arguments xdata, ydata, must be COLUMN vectors.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

sx = size(xdata);                   sy = size(ydata);
if sx(1) == 1 & sy(1) == 1 
    disp('Detected row vectors as input. Converting to column vectors.');
    xdata = xdata';                 ydata = ydata';
end

k = ko;        

% Call fminsearch with a random starting point.
start_point = rand(1);
model = @fitfun;
estimates = fminsearch(model, start_point);
% fitfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for A * exp(-lambda * xdata) - ydata, 
% and the FittedCurve. FMINSEARCH only needs sse, but we want to 
% plot the FittedCurve at the end.
    function [sse, FittedCurve] = fitfun(params)
        A = params;
        FittedCurve = A .* exp(-k * xdata);
        ErrorVector = FittedCurve - ydata;
        sse = sum(ErrorVector .^ 2);
    end

[sse, FittedCurve] = model(estimates);

% Calculate Correlation Coefficient
Fit = estimates .* exp(-k * xdata);
Compare = [ydata Fit];
[r p] = corrcoef(Compare);
R = r(1,2);     p_value = p(1,2);
% --------------------------------

% ---------- Plot -----------
% plot(xdata, ydata, 'xb','MarkerSize',4)
% hold on;
% plot(xdata, FittedCurve, 'r');
% xlabel('x');        ylabel('f(estimates,x)');
% title(['Fitting to function Ae^{-kx} for fixed k=' num2str(k)]);
% legend('data', ['y = Ae^{-kx}']);
% xrange = xlim; xpos = xrange(2);
% yrange = ylim; ypos = 3*yrange(2)/4;
% text(xpos,ypos,['k=' num2str(k) ' ; A=' num2str(estimates)...
%     ' R=' num2str(R) ', p=' num2str(p_value)], ...
%     'HorizontalAlignment','right');
% hold off

% lambda  = log(2)/k;     % half-life
% disp(['lambda = ' num2str(lambda)]);

end