function Rsq = CoefDet(data,theor)
% Find Coefficient of Determination, R^2
% This works for any curve fit, linear or nonlinear.
% Data and theoretical curve vectors should be of the same length.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

SST = sum((data - mean(data)).^2);        % Total sum of squares for simulation data
SSR = sum((data - theor).^2);             % sum of square residuals (data vs theoretical predictions)

Rsq = 1 - SSR ./ SST;                     % Definition of R^2

% R = sqrt(Rsq);                            % Correlation Coefficient R
