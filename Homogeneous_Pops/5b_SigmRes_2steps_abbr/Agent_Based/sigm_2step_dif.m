function [dydt] = sigm_2step_dif(t,y)
% Sigmoidal Response system of DEs
% R <--> Rp              Rates: k_f, k_r
%     ^
%     |
%     2S                 S: Signal (reactant in forward reaction)
% Assume number of S molecules does not change.
% Agent-based interpretation of S population size.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

global k_f k_r S;

R = y(1);       Rp = y(2);
dydt = zeros(2,1);

dydt(1) = - k_f * R * S*(S-1) + k_r * Rp;
dydt(2) = + k_f * R * S*(S-1) - k_r * Rp;
