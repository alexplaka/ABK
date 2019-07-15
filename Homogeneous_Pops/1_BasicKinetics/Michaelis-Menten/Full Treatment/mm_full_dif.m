function [dydt] = mm_dif(t,y)

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

global k_f k_r k_cat;

A = y(1);       B = y(2);       E = y(3);         EA = y(4);
dydt = zeros(4,1);

dydt(1) = - k_f * E * A + k_r * EA;
dydt(2) = + k_cat * EA;
dydt(3) = - k_f * E * A + k_r * EA + k_cat * EA;
dydt(4) = + k_f * E * A - k_r * EA - k_cat * EA;
