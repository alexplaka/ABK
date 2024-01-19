%  --> A
% Simulating 0th order kinetics using ABK.
% Generating a heterogeneous population of A with 
% normally-distributed k values (where k is the rate
% constant of A agents with respect to some unspecified
% and unsimulated process -- for instance, a degradation
% reaction in a birth-death process). 

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear; clc;
rng(0);
global k_b;               % Zeroth order MICROSCOPIC kinetic constant [units: N_a / sec];

k_b = 1;

agents = 100;          % Max number of A molecules that can be produced 
Ao = 0;                 % Initial number of A molecules
A(1) = Ao;              % Initialize time-dependent sum of A molecules

% Note: I don't need an array for tracking the status of individual A agents
% because that is not necessary since A is produced. That is, there won't be 
% a loop going through all A agents

time(1) = 0;
t = 1;                  % Time Counter variable
dt = 0.01;              % Fixed time step increment

Av = zeros(1,agents);

% Rate constant for each new agent A, for some subsequent reaction...
k = zeros(1,agents);        % generating heterogeneity

while A(t) < agents               

    P = k_b * dt;         % Does NOT depend on A
    
    if rand < P       
        A(t+1) = A(t) + 1;
        
        x = ceil(rand*agents);
        while Av(x) == 1        % ensure this agent hasn't been born already
            x = ceil(rand*agents);
        end
        Av(x) = 1;
        
        % create normally distributed k values with mean=0.6 and std=0.1
        k(x) = 0.6 + 0.1 * randn;  
        
    else
        A(t+1) = A(t);
    end
    time(t+1) = time(t) + dt;
    t = t + 1;
end

%% Plot time course
figure('Name','0th Order Rx Time course','NumberTitle','off'); 
scatter(time,A,3,'.r');                hold on;
xlabel('time');                 ylabel('A');           

% Solve ODE for 0th order kinetics - Plot Time Course
tmax=time(t-1);
[t_sol,y_sol] = ode45(@o0_dif,[0:tmax/100:tmax],Ao);
scatter(t_sol,y_sol(:,1),'.b');       
legend('ABK','DE');
