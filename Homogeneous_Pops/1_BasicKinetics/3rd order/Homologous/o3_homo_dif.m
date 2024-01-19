function [dydt] = o3_homo2_dif(t,y)
% Solve ODE for 3A --> D

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

global k dif_canon;               % MICROscopic kinetic constant
A = y(1);               D = y(2);

dydt = zeros(2,1);

if exist('dif_canon') == 0
    dif_canon = 0;                          % Default: Microscopic form of DE
end

if dif_canon == 0                               % Microscopic interpretation
    dydt(1) = - 3 * k * A * (A-1) * (A-2);      % This is consistent with the definition of
    dydt(2) = +     k * A * (A-1) * (A-2);      % the rx rate r = - 1/a * dNa/dt
elseif dif_canon == 1                           % Canonical form of DE
    dydt(1) = - 3 * k * A^3;
    dydt(2) = +     k * A^3;
end

% Keep in mind that for the reaction 3A --> D
% Note that dD/dt =  - 1/3 * dA/dt

% Also notice that k is the microscopic rate constant, and therefore
% A represents the number of molecules. Hence, the differential
% equation contains the term k * A * (A-1) * (A-2).  -- For the microscopic version --