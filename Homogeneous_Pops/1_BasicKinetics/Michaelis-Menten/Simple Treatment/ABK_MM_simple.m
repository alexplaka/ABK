% Simulating Michaelis-Menten kinetics - Simple Implementation
% Overall rx: A --> B
% Actual rxs: E + A <--> EA --> E + B       (Not explicitly modeled)
% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear;       tic;
rng(0);

global k_cat Km E_tot;

maxTime = 500;              % Max Simulation time (sec)
dt = 1/100;                   % Constant (fixed) time step increment (sec)     
t_steps = maxTime / dt;
time = zeros(1,t_steps);

k_cat = 0.05;
Km = 10;                   % MM constant: # of molecules required for half-maximal speed
E_tot = 5;                 % Total number of Enzyme molecules

agents = 20;
A = ones(2,agents);      
Sa = zeros(1,t_steps);           Sa(1) = agents;        % ** Initial condition **
v(1) = k_cat*E_tot*Sa(1) / (Sa(1)+Km); 
P = zeros(1,t_steps);
t = 2;

while Sa(t-1) >= 1 && t <= t_steps
    for i = 1:size(A,2)
        if A(1,i) == 1                                  % if agent is still alive
            P(t-1) = k_cat*E_tot / (Sa(t-1)+Km) * dt;   % P from MM kinetics rate law
            if rand < P(t-1)                            % **check probability condition**
                A(2,i) = 0;                             % A agent dies (converted to B)
            end
        end
    end
    time(t) = time(t-1) + dt;      
    Sa(t) = sum(A(2,:));
    v(t) = k_cat*E_tot*Sa(t) / (Sa(t)+Km);              % Classic MM kinetics; -dA/dt = dB/dt

    A(1,:) = A(2,:);
    t = t + 1;
end

% Remove unnecessary terminal 0's from arrays
Sa = Sa(Sa~=0);               P = P(1:size(Sa,2));          v = v(1:size(Sa,2));
time = time(1:size(Sa,2));    
Sb = agents - Sa;

%% Solve differential equation
finaltime = time(end);
[t_sol, y_sol] = ode45(@mm_sim_dif,0:maxTime/100:maxTime,[Sa(1) ; Sb(1)]);

%% Plot Process Time Course
figure;                                                             hold on;
scatter(time,Sa,'.m');                      scatter(time,Sb,'.g');
scatter(t_sol,y_sol(:,1),3,'om');           scatter(t_sol,y_sol(:,2),3,'og');
xlabel('t');                                ylabel('agents');        
legend('ABM A','ABM B','DE A','DE B');                              hold off;

%% plot Rx Speed/Rate vs A (hyperbolic)
figure;
plot(Sa,v,'r','LineWidth',3);                       
xlabel('N_A');    ylabel('-dA/dt');
axis([0 agents 0 k_cat*E_tot]);
%% plot P vs A
% plot(Sa,P,'r','LineWidth',3);                       
% xlabel('N_A');    ylabel('P_{dif,MM}');
% axis([0 agents 0 max(P)]);

%% Finish
clear t w i;
toc

% Note:
% - The simplicity of this implementation makes it ideal for
% simulating reactions obeying Michaelis-Menten kinetics,
% esp. when the population of A is large compared to E_tot.