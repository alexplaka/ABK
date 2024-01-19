function [dydt] = o12_mixed_dif(t,y)
%  A     --> C              Rate constant k1  (1st order)
%  A + B --> D              Rate constant k2  (2nd order)

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

global k1 k2;
A = y(1);       B = y(2);       
C = y(3);       D = y(4);           % or X, Y

dydt = zeros(4,1);

dydt(1) = -k1*A -k2*A*B;
dydt(2) = -k2*A*B;
dydt(3) = +k1*A;
dydt(4) = +k2*A*B;