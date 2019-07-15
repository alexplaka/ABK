function [dydt] = o23_mixed2_dif(t,y)
%  A + B     --> X              Rate constant k1  (2nd order)
%  A + C + D --> Y              Rate constant k2  (3rd order)

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

global k1 k2;
A = y(1);       B = y(2);       C = y(3);       D = y(4);       
X = y(5);       Y = y(6);

dydt = zeros(6,1);

dydt(1) = -k1*A*B -k2*A*C*D;
dydt(2) = -k1*A*B;
dydt(3) = -k2*A*C*D;
dydt(4) = -k2*A*C*D;
dydt(5) = +k1*A*B;
dydt(6) = +k2*A*C*D;