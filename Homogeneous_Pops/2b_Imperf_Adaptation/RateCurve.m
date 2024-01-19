% Note that changing the value of S (Signal) does NOT change the steady-state 
% value of R (Response). This agrees with the theoretically determined R_ss 
% value, which does not depend on S.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear;
% Initialize
agents = 1000;

N_avo = 6.02e23;        % Avogadro's number
V = 1e-21;              % Volume in L

k1 = 1;                  % Synthesis of R; 1st order [units: 1/sec] or 0th order
k2_m = 10;                       % MOLAR rate of degradation of R; 1st order [units: 1/(M*sec)]
k2 = k2_m / (N_avo * V);       % MICROscopic rate of degradation of R; 1st order [units: 1/sec]
k3 = 1;                      % Synthesis of X; 1st order [units: 1/sec]
k4 = 0.1;                        % Degradation of X; 1st order [units: 1/sec]
K = 50;             % K_50 for X inhibiting production of R (# of X to produce half-maximal effect) 
nH = 1;             % Hill coefficient for X inhibiting the production of R

S = 20;

% Construct Rate curve
R_array = 0:agents;
% Rate curve must be done for a particular value of X. 
% This is an arbitrarily chosen value for X, yet one 
% reflecting its overall dependence on S.
X = S;

dRdt_degr = k2 * R_array;               
dRdt_synt = k1 * S * K^nH/(K^nH + X^nH) * ones(1,size(R_array,2));

[R_ss_int, ssr] = intersections(R_array,dRdt_degr,R_array,dRdt_synt);
% R_ss_int: Steady-state value of R derived from intersections in Rate curve
% ssr: Steady-state rate of forward and back reactions
figure('Name','Rate Curve','NumberTitle','off');
plot(R_array,dRdt_degr,'r',R_array,dRdt_synt,'-k',R_ss_int,ssr,'ob');
xlabel('R');        ylabel('dR/dt');    
legend('degradation','synthesis','R_{ss}','Location','Best');
clear R_array;
