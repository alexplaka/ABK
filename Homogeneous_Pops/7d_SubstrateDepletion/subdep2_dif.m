function [dydt] = subdep2_dif(t,y)
%          Ep  <--->  E          Rates: k_f, k_r (Reverse rx activated by R)
%          |     ^               Michaelis-Menten constants: Km_f, Km_r
%           \     \              k_s: synthesis of X
%    k_s     \     |             k_b1,2: synthesis of R
%  S ---> X  --->  R  ------>    k_d: degradation of R                
%           k_b1,2      k_d               ** S: Signal, R: Response **  

% Simulating substrate depletion oscillator process: 
% k_s: X synthesis, 0th order process, UPregulated by S
% k_b1: X --> R, 1st order process wrt X
% k_b2: X --> R, 2nd order process wrt X, Ep (Ep is not consummed)
% k_d: R degradation, 1st order process wrt R
% k_f: E synthesis, MM process
% k_r: Ep synthesis, MM process wrt R (but R is not consummed)
% Km_f, Km_r: MM constants for forward and reverse rxs

% Note: Assume number of S molecules/agents is NOT changing (eg, S is DNA molecules)

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

global agents k_s k_b1 k_b2 k_d k_f k_r Km_f Km_r S;

R = y(1);         Ep = y(2);             X = y(3);

dydt = zeros(3,1);

dydt(1) = + k_b1 * X + k_b2 * Ep * X - k_d * R;
dydt(2) = + k_r * (agents - Ep) * R / (Km_r + (agents - Ep)) - k_f * Ep / (Km_f + Ep);
dydt(3) = + k_s * S - k_b1 * X - k_b2 * Ep * X;