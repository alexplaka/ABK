function [dydt] = mm_sim_dif(t,y)
% Simulating Michaelis-Menten kinetics - Simple Implementation
% Overall rx: A --> B

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

global k_cat Km E_tot;

A = y(1);       B = y(2);
dydt = zeros(2,1);

dydt(1) = - k_cat * E_tot * A / (A + Km);
dydt(2) = + k_cat * E_tot * A / (A + Km);

