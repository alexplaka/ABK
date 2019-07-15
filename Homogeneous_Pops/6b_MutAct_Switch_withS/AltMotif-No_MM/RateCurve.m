function [R_ss_int, Ep_ss] = RateCurve()

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

global agents k_b k_d k_s k_f k_r S;

% agents = 1000;
% k_b = 0.01;                        % 1st order R synthesis rate
% k_d = 0.02;                        % 1st order degradation rate (units: 1/sec)
% k_s = 0.1;                           % 1st order R synthesis rate (units: 1/sec)
% k_f = 2;                         % basal forward rate (1st order)
% k_r = 0.01;                         % basal reverse rate (1st order)
% S = 20;                            % Assume number of S molecules/agents is NOT changing

%%
syms E Ep R positive;
dEp_sym = + k_r * E * R - k_f * Ep ;
temp = subs(dEp_sym,E,agents-Ep);
Ep_sym = solve(temp==0,Ep);                 % result has 2 entries, +/-
% Denominator contains R - k_f/k_r
%% Construct Rate curve
R_array = 0:50;              
R_array(k_f/k_r + 1) = [];      % Remove R=k_f/k_r entry (avoid div by 0)
Ep_array = double(subs(Ep_sym(1),R,R_array));

dRdt_degr = k_d .* R_array;
dRdt_synt = k_b .* Ep_array + k_s * S; 

%%
[R_ss_int, ssr] = intersections(R_array,dRdt_degr,R_array,dRdt_synt);
% R_ss_int: Steady-state value of R derived from intersections in Rate curve
% ssr: Steady-state rate of forward and back reactions

figure('Name','Rate Curve','NumberTitle','off');
plot(R_array,dRdt_degr,'r',R_array,dRdt_synt,'-k',R_ss_int,ssr,'ob');
xlabel('R');        ylabel('dR/dt');    
legend('degradation','synthesis','R_{ss}','Location','Best');
clear R_array;

Ep_ss = double(subs(Ep_sym(1),R, R_ss_int));