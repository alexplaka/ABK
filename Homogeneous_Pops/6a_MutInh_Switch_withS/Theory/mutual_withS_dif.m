function [dydt] = mutual_withS_dif(t,y)

%    k_da           k_db
%    <--  A _  _ B  -->
%         ^  \/  ^
%    k_sa |  /\  | k_sb
%         | / ~\ |
%         |~  | ~|
%             |
%             S
% A and B are synthesized (0th order constants k_sa, k_sb respectively)
% A and B are degraded (k_da, k_db)
% A and B inhibit each other's synthesis (MUTUAL inhibition, alpha < 1):
% B influences the rate of the synthesis of A
% A influences the rate of the synthesis of B
% S inhibits the inhibition of B by A (assume complete repression)
% Assume S does not change.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

global k_da k_db;
global k_a alpha_a K_b n_b;
global k_b alpha_b K_a n_a;
global K_s n_s S;

A = y(1);       B = y(2);
dydt = zeros(2,1);

k_sa = k_a * (K_b^n_b + alpha_b * B^n_b) / (K_b^n_b + B^n_b);
k_sb = k_b * (K_a^n_a + alpha_a * (A * 1/(1+(S/K_s)^n_s))^n_a) / (K_a^n_a + (A * 1/(1+(S/K_s)^n_s))^n_a);

dydt(1) = + k_sa - k_da * A;
dydt(2) = + k_sb - k_db * B;
