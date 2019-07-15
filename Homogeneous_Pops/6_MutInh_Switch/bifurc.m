function [] = bifurc(A_i,B_i)                          
% Input: Initial population sizes of A and B

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

global agents;
global k_da k_db;
global k_a alpha_a K_b n_b;
global k_b alpha_b K_a n_a;

totalTime = 5000;                    % Simulation time (sec)

A_ini = 0:10:agents;
B_ini = B_i;

figure('Name','Deterministic Time Course','NumberTitle','off');       hold on;

for i=1:size(A_ini,2)
    [t_sol, y_sol] = ode45(@mutual_dif,0:totalTime/500:totalTime,[A_ini(i); B_ini]);
    scatter(t_sol,y_sol(:,1),3,'.b');
    A_fin(i) = y_sol(end,1);
end

xlabel('t (sec)');                   ylabel('N_A');

% Find at what value of A_ini bifurcation occurs
Rd = abs(diff(A_fin));
bf = find(Rd > agents/12);
if isempty(bf)==1,     disp('No bifurcations predicted');     return;     end
fprintf(['Bifurcation occurs for:\t' num2str(A_ini(bf(1))) ' < A_ini < ' num2str(A_ini(bf(end)+1)) '\n']);
sink(1) = mean(A_fin(1:bf(1)));
sink(2) = mean(A_fin(bf(end)+1:end));
fprintf(['Stable A_ss values: ' num2str(sink(1)) ' <-- A_fin --> ' num2str(sink(2)) '\n']);

% Looking more closely to find smaller range of R_ini leading to bifurcation
A_ini2 = A_ini(bf(1)):A_ini(bf(end)+1);

for j=1:size(A_ini2,2)
    [t_sol, y_sol] = ode45(@mutual_dif,0:totalTime/500:totalTime,[A_ini2(j); B_ini]);
    scatter(t_sol,y_sol(:,1),3,'.c');
    A_fin2(j) = y_sol(end,1);
end
Rd2 = abs(diff(A_fin2));
bf2 = find(Rd2 > agents/25);
disp('Looking more closely...'); 
fprintf(['Bifurcation occurs for:\t' num2str(A_ini2(bf2(1))) ...
    ' < A_ini < ' num2str(A_ini2(bf2(end)+1)) '\n\n']);

% plot trajectory for specific initial value of A
[t_sol, y_sol] = ode45(@mutual_dif,0:totalTime/500:totalTime,[A_i; B_ini]);
scatter(t_sol,y_sol(:,1),3,'.r');                                     


