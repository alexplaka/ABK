function [dydt] = negFb3_dif2(t,y)
%          R <---->  Rp             Rates: k_f, k_r (Forward rx activated by Yp)
%               ^     |             Michaelis-Menten constants: Km_f, Km_r
%              /      |
%             |       |
%    Y <---> Yp       |             Rates: k_f1, k_r1 (Forward rx activated by X)
%         ^           |             Michaelis-Menten constants: Km_f1, Km_r1
%        /            |
%        |            \
%   -->  X  ---------------->                      Rp: Response
%   k_b  ^       k_d1/k_d2
%        | 
%        | k_s
%        S                             S: Signal

% Simulating negative feedback process (3-component system: X, Y/Yp, R/Rp) 
% k_b: X synthesis, 0th order process
% k_s: X synthesis, 1st order process wrt S (but S is not consummed)
% k_s: ALTERNATIVE IMPLEMENTATION: X synthesis, 0th order process, UPregulated by S
% k_d1: X degradation, 1st order process wrt X
% k_d2: X degradation, 2nd order process wrt X, Rp (but Rp is not consummed)
% k_f: Rp synthesis, MM process wrt Yp (but Yp is not consummed)
% k_r: R synthesis, MM process 
% Km_f, Km_r: MM constants for forward and reverse rxs R <==> Rp
% k_f1: Yp synthesis, MM process wrt X (but X is not consummed)
% k_r1: Y synthesis, MM process 
% Km_f1, Km_r1: MM constants for forward and reverse rxs Y <==> Yp

% Note: Assume number of S molecules/agents is NOT changing (eg, S is DNA molecules)

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

global  k_s k_b k_d1 k_d2 k_f k_r Km_f Km_r k_f1 k_r1 Km_f1 Km_r1 S;

% Model 5 species (X, Y, Yp, R, Rp)
X = y(1);         
Y = y(2);               Yp = y(3);            
R = y(4);               Rp = y(5);
dydt = zeros(5,1);

dydt(1) = + k_b + k_s * S - k_d1 * X - k_d2 * X * Rp;
dydt(2) = - k_f1 * Y * X / (Km_f1 + Y) + k_r1 * Yp / (Km_r1 + Yp);
dydt(3) = + k_f1 * Y * X / (Km_f1 + Y) - k_r1 * Yp / (Km_r1 + Yp);
dydt(4) = - k_f  * R * Yp / (Km_f  + R) + k_r * Rp / (Km_r + Rp);
dydt(5) = + k_f  * R * Yp / (Km_f  + R) - k_r * Rp / (Km_r + Rp);




