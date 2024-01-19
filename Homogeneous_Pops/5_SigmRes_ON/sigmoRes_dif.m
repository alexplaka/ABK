function [dydt] = sigmoRes_dif(t,y)
% Sigmoidal Response system of DEs: 
% (Michaelis-Menten kinetics for forward and reverse rxs)
% R <--> Rp             Rates: k_f, k_r
%     ^                 Michaelis-Menten constants: Km_f, Km_r
%     |
%     S                 S: Signal (reactant in forward reaction)
% Assume number of S molecules does not change.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

global k_f k_r Km_f Km_r S;

R = y(1);       Rp = y(2);
dydt = zeros(2,1);

dydt(1) = - k_f * R * S / (Km_f + R)  +  k_r * Rp / (Km_r + Rp);
dydt(2) = + k_f * R * S / (Km_f + R)  -  k_r * Rp / (Km_r + Rp);
