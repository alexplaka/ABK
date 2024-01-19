function [dydt] = bruss_can_dif(t,y)
% The Brusselator rx system is:
%   --> X           (birth of X - 0th order; rate constant: b)
% X -->             (death or degradadation of X - 1st order wrt X; rate constant: d)
% X --> Y           (1st order wrt X; rate constant: c)
% 2X + Y --> 3X     (3rd order: 2nd wrt X, 1st wrt Y; rate constant: a)

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

global a b c d;

X = y(1);                   Y = y(2);

dydt = zeros(2,1);

% ** Canonical formulation of DEs **
dydt(1) = b + a * X^2 * Y - (c + d) * X;        
dydt(2) = c * X - a * X^2 * Y;         
