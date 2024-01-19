function [dydt] = actinh_dif(t,y)
%   Ep <--->  E            Rates: k_f, k_r (Reverse rx activated by R)
%    |  ^                  Michaelis-Menten constants: Km_f, Km_r
%    |   \   --> X -->     k_bx and k_dx: synthesis and degradation of X
%    \   \  /    \
%    -->  R  -------->                      R: Response
%   k_b   ^      k_d
%         | k_s
%         S                                 S: Signal

% Simulating mutual activation switch process for R: 
% k_b: R synthesis, 0th order process (Ep is not consummed)
% k_d1: R degradation, 1st order process wrt R
% k_d2: R degradation, 2nd order process wrt R, X
% k_bx: X synthesis, 0th order, UPregulated by R
% k_dx: X degradation, 1st order wrt X
% k_s: R synthesis, 1st order process wrt S
% k_f: E synthesis, MM process
% k_r: Ep synthesis, MM process wrt R (but R is not consummed)
% Km_f, Km_r: MM constants for forward and reverse rxs

% Note: Assume number of S molecules/agents is NOT changing (eg, S is DNA molecules)

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

global k_b k_d1 k_d2 k_s k_bx k_dx k_f k_r Km_f Km_r S;
R = y(1);         Ep = y(2);          E = y(3);         X = y(4);
dydt = zeros(4,1);

dydt(1) = + k_b * Ep + k_s * S - k_d1 * R - k_d2 * R * X;
dydt(2) = + k_r * E * R / (Km_r + E) - k_f * Ep / (Km_f + Ep);
dydt(3) = - k_r * E * R / (Km_r + E) + k_f * Ep / (Km_f + Ep);
dydt(4) = + k_bx * R - k_dx * X;