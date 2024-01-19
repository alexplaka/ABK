function [dydt] = adapt_dif(t,y)
%    k3     k4
%    -->  X -->                        
%   /     /
%   S    /                              S: Signal
%   \   /
%    -->  R -->                         R: Response
%    k1     k2

% Simulating birth-death process for R and X, and X inhibits the production of R: 
% k1: synthesis of R, 0th order process (upregulated by S; inhibited by X)
% k2: degradation of R, 1st order process (R)
% k3: synthesis of X, 0th order process (upregulated by S)
% k4: degradation of X, 1st order process (X)

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

global k1 k2 k3 k4 S K nH;

R = y(1);           X = y(2);
dydt = zeros(2,1);

dydt(1) = k1 * S * K^nH/(K^nH + X^nH) - k2 * R;
dydt(2) = k3 * S - k4 * X;


