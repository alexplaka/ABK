function [R_ss_int, E_ss] = RateCurve()

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

global maxRagents agents k_b k_d k_f k_r Km_f Km_r S;

%
syms E Ep R positive;
dEp_sym = + k_r * E * R / (Km_r + E) - k_f * Ep / (Km_f + Ep);
temp = subs(dEp_sym,E,agents-Ep);
Ep_sym = solve(temp==0,Ep);                 % result has 2 entries, -/+
% Denominator contains R - k_f/k_r

% Construct Rate curve
R_array = 0:maxRagents;     
k = k_f / k_r;
if k < maxRagents 
    R_array(k_f/k_r + 1) = [];      % Remove R=k_f/k_r entry (avoid div by 0)
end
Ep_array = double(subs(Ep_sym(2),R,R_array));
E_array = agents - Ep_array;

dRdt_degr = k_d .* R_array * S;
dRdt_synt = k_b .* E_array; 

%
[R_ss_int, ssr] = intersections(R_array,dRdt_degr,R_array,dRdt_synt);
% R_ss_int: Steady-state value of R derived from intersections in Rate curve
% ssr: Steady-state rate

% figure('Name','Rate Curve','NumberTitle','off');
% plot(R_array,dRdt_degr,'r',R_array,dRdt_synt,'-k',R_ss_int,ssr,'ob');
% xlabel('R');        ylabel('dR/dt');    
% legend('degradation','synthesis','R_{ss}','Location','Best');

Ep_ss = double(subs(Ep_sym(1),R, R_ss_int));
E_ss = agents - Ep_ss;