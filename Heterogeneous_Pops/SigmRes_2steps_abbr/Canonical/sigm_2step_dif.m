function [dydt] = sigm_2step_dif(t,y)
% Sigmoidal Response system of DEs
% R <--> Rp              Rates: k_f, k_r
%     ^
%     |
%     2S                 S: Signal (reactant in forward reaction)
% Assume number of S molecules does not change.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

global k_f k_r_mean S;

R = y(1);       Rp = y(2);
dydt = zeros(2,1);

dydt(1) = - k_f * R * S^2 + k_r_mean * Rp;
dydt(2) = + k_f * R * S^2 - k_r_mean * Rp;
