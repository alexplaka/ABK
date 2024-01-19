function [dydt] = negFb_delaydif(t,y,Z)
%    R  <--->  Rp            Rates: k_f, k_r (Forward rx activated by X)
%          ^    |            Michaelis-Menten constants: Km_f, Km_r
%         /     |
%         |     \
%    -->  X    -->                      Rp: Response
%   k_b   ^    k_d1/k_d2
%         | k_s
%         S                             S: Signal

% Simulating negative feedback process: 
% k_b: X synthesis, 0th order process
% k_s: X synthesis, 1st order process wrt S (but S is not consummed)
% k_s: ALTERNATIVE IMPLEMENTATION: X synthesis, 0th order process, UPregulated by S
% k_d1: X degradation, 1st order process wrt X
% k_d2: X degradation, 2nd order process wrt X, Rp (but Rp is not consummed)
% k_f: Rp synthesis, MM process wrt X (but X is not consummed)
% k_r: R synthesis, MM process 
% Km_f, Km_r: MM constants for forward and reverse rxs

% Note: Assume number of S molecules/agents is NOT changing (eg, S is DNA molecules)

% Using delay in feedback reaction.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

global agents k_s k_b k_d1 k_d2 k_f k_r Km_f Km_r S;

% Implementation: model two species (X, Rp)             PREFERRED

X = y(1);         Rp = y(2);       Rp_lag = Z(:,1);

dydt = zeros(2,1);

dydt(1) = + k_b + k_s * S - k_d1 * X - k_d2 * X * Rp_lag(2);
dydt(2) = - k_r * Rp / (Km_r + Rp) + k_f * (agents-Rp) * X / (Km_f + (agents-Rp));
