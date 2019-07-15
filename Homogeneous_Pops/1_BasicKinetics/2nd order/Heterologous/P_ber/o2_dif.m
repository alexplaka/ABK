function [dydt] = o2_dif(t,y)

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

global k
A = y(1);       B = y(2);       C = y(3);

dydt = zeros(3,1);

dydt(1) = -k*A*B;
dydt(2) = -k*A*B;
dydt(3) = +k*A*B;