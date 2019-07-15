function [dydt] = bd_dif(t,y)

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

global k_b k_d;

A = y(1);
dydt = zeros(1,1);

dydt(1) = k_b - k_d * A;
