function [dydt] = o3_het_dif(t,y)

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

global k
A = y(1);       B = y(2);       C = y(3);       D = y(4);

dydt = zeros(4,1);

dydt(1) = -k*A*B*C;
dydt(2) = -k*A*B*C;
dydt(3) = -k*A*B*C;
dydt(4) = +k*A*B*C;