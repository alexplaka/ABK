function [dydt] = hyperAdapt_dif(t,y)
% Hyperbolic Response Motif --> Adaptation: system of DEs
% 
%  S
% --> X -->              Rates: k_b*S, k_d
%     |
%     / 
% R <---> Rp             Rates: k_f, k_r
%     /
%     |
%     S                 S: Signal (reactant in forward reaction)

% Assume number of S molecules does not change.

% k_f: synthesis of Rp, 1st order process wrt R (upregulated by S), or 2nd order wrt R,S
% k_r: synthesis of R, 2nd order process (X, Rp), but X is not used up.
% k_b: synthesis of X, 0th order process (upregulated by S)
% k_d: degradation of X, 1st order process (X)

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

global k_b k_d k_f k_r S;

R = y(1);       Rp = y(2);      X = y(3);
dydt = zeros(3,1);

dydt(1) = - k_f * R * S + k_r * Rp * X;
dydt(2) = + k_f * R * S - k_r * Rp * X;
dydt(3) = + k_b * S     - k_d * X;
