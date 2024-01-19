function [dydt] = hadapt_dif(t,y)
% Hyperbolic Response - Adaptation system of DEs

%             |                        Rate: k_sp
%             \
% --> R <---> Rp -->                   Rates: k_sr, k_f, k_r, k_dp
%        /      /
%        |     /                       S: Signal (reactant in R forward reaction)
%        S    /             
%        |    |
%        \    |                        S: Signal (upregulates synthesis of X)
%        -->  X  -->               
%        k_sx    k_dx                  Synthesis and degradation rates of X

% k_sr:  0th order
% k_f:   2nd order wrt S, R
% k_r:   1st order wrt Rp
% k_sp:  0th order
% k_dp:  2nd order wrt Rp, X
% k_sx:  0th order (upregulated by S)
% k_dx:  1st order wrt X

% Assume number of S molecules does NOT change.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

global k_sr k_f k_r k_sp k_dp k_sx k_dx S;

R = y(1);       Rp = y(2);      X = y(3);
dydt = zeros(3,1);

dydt(1) = + k_sr - k_f * R * S + k_r * Rp;
dydt(2) = + k_sp + k_f * R * S - k_r * Rp - k_dp * Rp * X;
dydt(3) = + k_sx * S - k_dx * X;