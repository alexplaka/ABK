function [dydt] = KM_dif(t,y)
% H + S --> S + S                       Rate constant k_c  (2nd order)
%     S --> D                           Rate constant k_d  (1st order)

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

global k_c k_d;

H = y(1);       S = y(2);       D = y(3);

dydt = zeros(3,1);

dydt(1) = - k_c * H * S;
dydt(2) = + k_c * H * S - k_d * S;
dydt(3) =               + k_d * S;
