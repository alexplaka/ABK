function [dydt] = mAct_dif(t,y)
%   Ep  <--->  E            Rates: k_f, k_r (Reverse rx activated by R)
%    |   ^                  
%    |    \
%    \    |
%    -->  R   -->                      R: Response
%   k_b   ^   k_d
%         | k_s
%         S                            S: Signal

% Simulating mutual activation switch process for R: 
% k_b: R synthesis, 1st order process wrt Ep (but Ep is not consummed)
% k_d: R degradation, 1st order process wrt R
% k_s: R synthesis, 1st order process wrt S
% k_f: E synthesis
% k_r: Ep synthesis (but R is not consummed)

% Note: Assume number of S molecules/agents is NOT changing (eg, S is DNA molecules)

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

global k_b k_d k_s k_f k_r S;
R = y(1);         Ep = y(2);          E = y(3);
dydt = zeros(3,1);

dydt(1) = + k_b * Ep + k_s * S - k_d * R;
dydt(2) = + k_r * E * R - k_f * Ep;
dydt(3) = - k_r * E * R + k_f * Ep;
