function [dydt] = predprey_dif(t,y)
% Differential equations for Volterra model.
%     R --> 2R              Density-dependent Birth of R; rate constant: a (1st order)
% F + R -->  F              Death of R; rate constant: b (2nd order)
% F     -->                 Death of F; rate constant: c (1st order)
% F + R --> 2F + R          Birth of F; rate constant: d (2nd order)

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

global a b c d K;

R = y(1);       F = y(2);
dydt = zeros(2,1);

dydt(1) = a*R*(1-R/K) - b*F*R;
dydt(2) = d*R*F - c*F;
