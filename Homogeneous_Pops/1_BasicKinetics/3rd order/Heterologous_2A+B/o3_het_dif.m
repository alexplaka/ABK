function [dydt] = o3_het_dif(t,y)
% Solve ODEs for 2A + B --> D

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3


global k;                       % MICROscopic kinetic constant
A = y(1);       B = y(2);       D = y(3);

dydt = zeros(3,1);

% dydt(1) = - 2 * k * A^2 * B;              % Canonical DE
% dydt(2) = -     k * A^2 * B;
% dydt(3) = +     k * A^2 * B;

dydt(1) = - 2 * k * A * (A-1) * B;          % Agent-Based DE
dydt(2) = -     k * A * (A-1) * B;
dydt(3) = +     k * A * (A-1) * B;


% Note that dB/dt =   1/2 * dA/dt < 0
% Note that dD/dt = - 1/2 * dA/dt > 0