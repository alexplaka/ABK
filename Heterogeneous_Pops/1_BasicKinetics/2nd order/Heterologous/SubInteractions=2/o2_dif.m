function [dydt] = o2_dif(t,y,k)

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

A = y(1);       B = y(2);       C = y(3);

dydt = zeros(3,1);

dydt(1) = -k * A * B;
dydt(2) = -k * A * B;
dydt(3) = +k * A * B;