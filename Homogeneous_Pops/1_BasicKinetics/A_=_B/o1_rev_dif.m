function [dydt] = o1_rev_dif(t,y)

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

global ko_f ko_r;

A = y(1);       B = y(2);
dydt = zeros(2,1);

dydt(1) = -ko_f * A + ko_r * B;
dydt(2) =  ko_f * A - ko_r * B;
