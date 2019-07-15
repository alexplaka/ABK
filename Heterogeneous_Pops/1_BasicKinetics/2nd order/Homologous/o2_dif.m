function [dydt] = o2_dif(t,y,k)

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

global DE_canon;

A = y(1);                               C = y(2);       

dydt = zeros(2,1);

if exist('DE_canon','var') == 0
    DE_canon = 1;                      % Default: Canonical form of DE
end

if DE_canon == 0
    dydt(1) = - 2 * k * A * (A-1);      % Agent-based form of DE
    dydt(2) = +     k * A * (A-1);      
elseif DE_canon == 1
    dydt(1) = - 2 * k * A^2;            % Canonical form of DE
    dydt(2) = +     k * A^2;      
end

% Keep in mind that for the reaction 2A --> C
% dC/dt = - 1/2 dA/dt

% Also notice that k is the microscopic rate constant, and therefore
% A represents the number of molecules. Hence, the differential
% equation contains the term k * A * (A-1).   -- For the microscopic version --