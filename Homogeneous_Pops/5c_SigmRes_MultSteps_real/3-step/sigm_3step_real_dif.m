function [dydt] = sigm_3step_real_dif(t,y)
% Sigmoidal Response system of DEs - 3-step process (more realistic)

% R <--> Rp <--> Rpp <--> Rppp           Rates: k_f, k_r
%     ^       ^        ^
%     |       |        |
%     S       S        S                 S: Signal (reactant in forward reactions)

% Assume number of S molecules does not change.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

global k_f k_r S;

R = y(1);       Rp = y(2);          Rpp = y(3);     Rppp = y(4);
dydt = zeros(4,1);

dydt(1) = - k_f * R * S + k_r * Rp;
dydt(2) = + k_f * R * S - k_r * Rp - k_f * Rp * S + k_r * Rpp;
dydt(3) = + k_f * Rp * S - k_r * Rpp - k_f * Rpp * S + k_r * Rppp;
dydt(4) = + k_f * Rpp * S - k_r * Rppp;
