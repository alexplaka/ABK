function [dydt] = hyperRes_dif(t,y)
% Hyperbolic Response system of DEs
% R <--> Rp             Rates: k_f, k_r
%     ^
%     |
%     S                 S: Signal (promotes forward reaction)
% Assume number of S molecules does not change.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

global k_f k_r S;

R = y(1);       Rp = y(2);
dydt = zeros(2,1);

dydt(1) = - k_f * R * S + k_r * Rp;
dydt(2) = + k_f * R * S - k_r * Rp;
