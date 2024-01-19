function [dydt] = bruss_ab_dif(t,y)
% The Brusselator rx system is:
%   --> X           (birth of X - 0th order; rate constant: b)
% X -->             (death or degradadation of X - 1st order wrt X; rate constant: d)
% X --> Y           (1st order wrt X; rate constant: c)

% Replace the 3rd order reaction
% 2X + Y --> 3X     (3rd order: 2nd wrt X, 1st wrt Y; rate constant: a)
% with: 
% 2X <--> Z                     k_f/k_r = K_eq
% Z + Y --> Z + X               2nd order rate constant a 

% Agent-Based interpretation

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

global a b c d k_f k_r;

X = y(1);                   Y = y(2);               Z = y(3);

dydt = zeros(3,1);

% ** Agent-Based formulation of DEs **
dydt(1) = b - 2 * k_f * X * (X-1) + 2 * k_r * Z + a * Z * Y  - (c + d) * X;        
dydt(2) = c * X - a * Z * Y;
dydt(3) = k_f * X * (X-1) - k_r * Z;


% Notice that when the last reaction is in equilibrium (steady-state), then
% Z* = (k_f / k_r) * X * (X-1)
% Substituting that into the first two DEs and simplifying, we get:
% 
% dydt(1) = b + a * (k_f / k_r) * X * (X-1) * Y - (c + d) * X;        
% dydt(2) = c * X - a * (k_f / k_r) * X * (X-1) * Y;
% 
% In other words, the system reduces to the unmodified agent-based Brusselator 
% with the rate constant for the trimolecular rx being a * (k_f / k_r).
