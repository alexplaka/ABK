function [] = bifurc(R_i,Ep_i)                          
% Input: Initial values of R and Ep

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

global agents k_b k_d k_s k_f k_r Km_f Km_r S;

% N_avo = 6.02e23;        % Avogadro's number
% V = 1e-21;              % Volume in L

totalTime = 4500;                    % Simulation time (sec)

% agents = 1000;
% k_b = 0.01;                        % 1st order R synthesis rate
% k_d = 0.02;                        % 1st order degradation rate (units: 1/sec)
% k_s = 0.1;                         % 1st order R synthesis rate (units: 1/sec)
% k_f = 2;                           % basal forward rate (1st order)
% k_r = 0.01;                        % basal reverse rate (1st order)
% Km_fm = 0.05;                      % MOLAR Michaelis-Menten constant for forward rx
% Km_rm = 0.5;                       % MOLAR Michaelis-Menten constant for reverse rx
% Km_f = Km_fm * N_avo * V;          % MICROSCOPIC Michaelis-Menten constant for forward rx
% Km_r = Km_rm * N_avo * V;          % MICROSCOPIC Michaelis-Menten constant for reverse rx
% S = 20;                            % Assume number of S molecules/agents is NOT changing

R_ini = 0:10:agents;
Ep_ini = Ep_i;
E_ini = agents - Ep_ini;

figure('Name','Deterministic Time Course','NumberTitle','off');       hold on;

for i=1:size(R_ini,2)
    [t_sol, y_sol] = ode45(@mAct_dif,0:totalTime/500:totalTime,[R_ini(i); Ep_ini; E_ini]);
    scatter(t_sol,y_sol(:,1),3,'.b');
    R_fin(i) = y_sol(end,1);
end

% plot trajectory for specific initial value of R
[t_sol, y_sol] = ode45(@mAct_dif,0:totalTime/500:totalTime,[R_i; Ep_ini; E_ini]);
scatter(t_sol,y_sol(:,1),3,'.r');                                     

xlabel('time (sec)');                   ylabel('# R');

% Find at what value of R_ini bifurcation occurs
Rd = abs(diff(R_fin));
bf = find(Rd > 1);
disp(['Bifurcation occurs for:    ' num2str(R_ini(bf)) ' < R_ini < ' num2str(R_ini(bf+1))]);
sink(1) = mean(R_fin(1:bf));
sink(2) = mean(R_fin(bf+1:end));
disp(['Stable R_ss values: ' num2str(sink(1)) ' <-- R_fin --> ' num2str(sink(2))]);

% Looking more closely to find smaller range of R_ini leading to bifurcation
R_ini2 = R_ini(bf):R_ini(bf+1);

for j=1:size(R_ini2,2)
    [t_sol, y_sol] = ode45(@mAct_dif,0:totalTime/500:totalTime,[R_ini2(j); Ep_ini; E_ini]);
    scatter(t_sol,y_sol(:,1),3,'.g');
    R_fin2(j) = y_sol(end,1);
end
Rd2 = abs(diff(R_fin2));
bf2 = find(Rd2 > 100);
disp('Looking more closely...'); 
disp(['Bifurcation occurs for:    ' num2str(R_ini2(bf2)) ' < R_ini < ' num2str(R_ini2(bf2+1))]);



