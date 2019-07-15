function [estimates, FittedCurve] = FitCurve_exp(xdata, ydata)
% Fits data to equation y = A e^(-kx)
% Parameters A and k are optimized to fit the curve and reported.

% Author: Alex Plakantonakis
% License: GNU GPLv3

% Input arguments xdata, ydata, should be COLUMN vectors.

sx = size(xdata);                   sy = size(ydata);
if sx(1) == 1 & sy(1) == 1 
    disp('Detected row vectors as input. Converting to column vectors.');
    xdata = xdata';                 ydata = ydata';
end

% Call fminsearch with a random starting point.
start_point = rand(1, 2);
model = @fitfun;
estimates = fminsearch(model, start_point);
% fitfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for A * exp(-lambda * xdata) - ydata, 
% and the FittedCurve. FMINSEARCH only needs sse, but we want to 
% plot the FittedCurve at the end.
    function [sse, FittedCurve] = fitfun(params)
        A = params(1);
        k = params(2);
        FittedCurve = A .* exp(-k * xdata);
        ErrorVector = FittedCurve - ydata;
        sse = sum(ErrorVector .^ 2);
    end

[sse, FittedCurve] = model(estimates);

% Calculate Correlation Coefficient
Fit = estimates(1) .* exp(-estimates(2) * xdata);
Compare = [ydata Fit];
[r p] = corrcoef(Compare);
R = r(1,2);     p_value = p(1,2);
% ---------------------------

% ---------- Plot -----------
% plot(xdata, ydata, '*');                            hold on;
% plot(xdata, FittedCurve, 'r')
 
% xlabel('x')
% ylabel('f(estimates,x)')
% title(['Fitting to function A exp^(-k * x)']);
% legend('data', ['y = A exp(-k * x)']);
% xrange = xlim; xpos = xrange(2);
% yrange = ylim; ypos = 3*yrange(2)/4;

% text(xpos,ypos,['A=' num2str(estimates(1)) ' ; k=' num2str(estimates(2))...
%     ' R=' num2str(R) ', p=' num2str(p_value)], ...
%     'HorizontalAlignment','right');
% hold off

end