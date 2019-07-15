function [dydt] = mixed_o2Rxs_dif(t,y)
%  A + B --> X              Rate constant k1  (2nd order)
%  A + C --> Y              Rate constant k2  (2nd order)

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

global k1 k2;
A = y(1);       B = y(2);       C = y(3);
X = y(4);       Y = y(5);

dydt = zeros(5,1);

dydt(1) = -k1*A*B -k2*A*C;
dydt(2) = -k1*A*B;
dydt(3) = -k2*A*C;
dydt(4) = +k1*A*B;
dydt(5) = +k2*A*C;