function [dydt] = mAct_dif(t,y)
%   E   <--->  Ep             Rates: k_r, k_f 
%    |     ^                  Michaelis-Menten constants: Km_f, Km_r
%    |    /                   Reverse rx E --> Ep activated by R
%    \    |   k_d
%    -->  R   -->                      R: Response
%   k_b        ^   
%              |
%              S                       S: Signal

% Simulating homeostatic process for R: 
% k_b: R synthesis, 1st order process wrt Ep (but Ep is not consummed)
% k_d: R degradation, 2nd order process wrt R, S
% k_f: E synthesis, MM process
% k_r: Ep synthesis, MM process wrt R (but R is not consummed)
% Km_f, Km_r: MM constants for forward and reverse rxs

% Note: Assume number of S molecules/agents is NOT changing (eg, S is DNA molecules)

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

global k_b k_d k_f k_r Km_f Km_r S;
R = y(1);         Ep = y(2);          E = y(3);
dydt = zeros(3,1);

dydt(1) = + k_b * E - k_d * R * S;
dydt(2) = + k_r * E * R / (Km_r + E) - k_f * Ep / (Km_f + Ep);
dydt(3) = - k_r * E * R / (Km_r + E) + k_f * Ep / (Km_f + Ep);
