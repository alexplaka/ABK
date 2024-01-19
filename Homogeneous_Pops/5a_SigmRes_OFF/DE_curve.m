% Sigmoidal Response system of DEs: 
% (Michaelis-Menten kinetics for forward and reverse rxs)
% R <--> Rp             Rates: k_f, k_r
%    ^                 Michaelis-Menten constants: Km_f, Km_r
%    |
%    S                 S: Signal (reactant in reverse reaction)
% Assume number of S molecules does not change.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

global k_f k_r Km_f Km_r S;

maxTime = 2500;                    % Maximum Simulation time (sec)

agents = 100;               % disp(['Agents = ' num2str(agents)]);

k_f = 1;                             % basal forward rate constant
k_r = 0.05;                              % basal reverse rate constant
Km_f = 5;                              % MICROSCOPIC MM constant for forward rx
Km_r = 30;                              % MICROSCOPIC MM constant for reverse rx

S = 20;                            % Number of S molecules is constant

R_o = floor(agents/2);       		Rp_o = agents - R_o;		

% Solve differential equation
[t_sol, y_sol] = ode45(@sigmoRes_dif,0:maxTime/100:maxTime,[R_o ; Rp_o]);
disp(['Diff eq terminal R value   = ' num2str(y_sol(end,1))]);
disp(['Diff eq terminal Rp value  = ' num2str(y_sol(end,2))]);

scatter(t_sol,y_sol(:,2),3,'.r');                               hold off;
axis([0 t_sol(end) 0 agents]);          xlabel('time');         ylabel('# Agents');         
