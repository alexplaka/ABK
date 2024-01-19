function [dydt] = o1f_simple_dif(t,y)

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

global ko alpha n K;

A = y(1);       B = y(2);
dydt = zeros(2,1);

if alpha == 1
    H = 0;                            % No Regulation
else
    H =  (B^n) / (K^n + B^n);         % Hill Function
end

k = ko + ko * (alpha - 1) * H;

dydt(1) = -k*A;
dydt(2) =  k*A;
