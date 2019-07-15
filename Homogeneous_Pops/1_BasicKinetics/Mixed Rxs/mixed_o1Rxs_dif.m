function [dydt] = mixed_o1Rxs_dif(t,y)
%  A     --> C              Rate constant k1  (1st order)
%  A     --> D              Rate constant k2  (1st order)

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

global k1 k2;
A = y(1);       C = y(2);       D = y(3);

dydt = zeros(3,1);

dydt(1) = -k1*A -k2*A;
dydt(2) = +k1*A;
dydt(3) = +k2*A;