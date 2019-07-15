function [dydt] = o1f_simple_dif(t,y)

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

global ko alpha n K;

A = y(1);       B = y(2);
dydt = zeros(2,1);

H =  (K^n + alpha * B^n) / (K^n + B^n);         % RSF
k = ko * H;

dydt(1) = -k*A;
dydt(2) =  k*A;
