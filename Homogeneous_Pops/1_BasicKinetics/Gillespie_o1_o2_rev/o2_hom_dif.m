function [dydt] = o2_hom_dif(t,y)

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

global k
A = y(1);            C = y(2);

dydt = zeros(2,1);

dydt(1) = -k*A*(A-1);
dydt(2) = +k*A*(A-1);