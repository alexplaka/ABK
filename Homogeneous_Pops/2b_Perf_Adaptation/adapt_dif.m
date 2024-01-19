function [dydt] = adapt_dif(t,y)
%    k3     k4
%    -->  X -->                        
%   /     |
%   S     |                             S: Signal
%   \     \
%    -->  R -->                         R: Response
%    k1     k2

% Simulating birth-death process for R and X, and their synthesis from S: 
% k1: synthesis of R, 0th order process (upregulated by S)
% k2: degradation of R, 2nd order process (X, R) but X is not consummed
% k3: synthesis of X, 0th order process (upregulated by S)
% k4: degradation of X, 1st order process (X)

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

global k1 k2 k3 k4 S;

R = y(1);           X = y(2);
dydt = zeros(2,1);

dydt(1) = k1 * S - k2 * X * R;
dydt(2) = k3 * S - k4 * X;


