function [dydt] = o2_homo_dif(t,y)
% Solving ODE for 2A --> C

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

A = y(1);                       C = y(2);
dydt = zeros(2,1);

global k dif_canon;

if exist('dif_canon','var') == 0
    dif_canon = 0;                      % Default: Microscopic (Agent-based) form of DE
end

if dif_canon == 0
    dydt(1) = - 2 * k * A * (A-1);      % Microscopic interpretation
    dydt(2) = +     k * A * (A-1);      
elseif dif_canon == 1
    dydt(1) = - 2 * k * A^2;            % Canonical form of DE
    dydt(2) = +     k * A^2;      
end


% Keep in mind that for the reaction 2A --> C
% dC/dt = - 1/2 dA/dt

% Also notice that k is the microscopic rate constant, and therefore
% A represents the number of molecules. Hence, the differential
% equation contains the term k * A * (A-1).   -- For the microscopic version --