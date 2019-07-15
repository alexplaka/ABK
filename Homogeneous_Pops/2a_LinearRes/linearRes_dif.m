function [dydt] = linearRes_dif(t,y)

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

global k_b k_d k_s S;

R = y(1);

dydt = zeros(1,1);

dydt(1) = k_b + k_s * S - k_d * R;


