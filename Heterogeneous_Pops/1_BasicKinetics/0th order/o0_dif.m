function [dydt] = o0_dif(t,y)

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

global k_b;

dydt = k_b;