%  --> A -->  
%  k_b   k_d
% Simulating birth-death process for A: 
% k_b: birth, 0th order process; k_d: death, 1st order process. 
% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear;       tic;
rng(0);

N_avo = 6.02e23;        % Avogadro's number
V = 1e-21;              % Volume in L
global k_b k_d;

k_bm = 0.001;                 % Molar "birth" rate; 0th order (units: M/sec)
k_b = k_bm * N_avo * V;       % Microscopic 0th order birth rate
k_d = 0.01;                   % "Death" rate; 1st order (units: 1/sec)

A_ss = k_b / k_d;             % Steady-state value
disp(['Predicted A_ss = ' num2str(A_ss)]);

agents = floor(2 * A_ss);

totalTime = 700;              % Simulation time (sec)
dt = 1/100;                   % Constant (fixed) time step increment (sec)

if exist('dt','var') == 1         
    t_steps = totalTime / dt;
else
    t_steps = 1200;             % empirically chosen value (for lambda = 0.8)
end
time = zeros(1,t_steps);

Sa = zeros(1,t_steps);              Sa(1) = 0;        % ** Initial condition **
tempA = zeros(1,agents);             
% Put "1" where agents are "alive", then randomize the array
for c=1:Sa(1),                  tempA(c)=1;             end
tempA = RandArray(tempA);                   % Randomize array
% Markov process, so only previous and current time steps needed --> 2 rows:    
A = [tempA ; tempA];                        % Initialize vector storing state of A agents
clear c tempA;

% P_b = k_b * dt;               % Probability of "birth" process (0th order), P_dif
P_b = 1 - exp(-k_b*dt);       % Probability of "birth" process (0th order), P_ber

% P_d = k_d * dt;               % Probability of "death" process (1st order), P_dif
P_d = 1 - exp(-k_d * dt);     % Probability of "death" process (1st order), P_ber


for t = 2:t_steps
     
    if rand < P_b                           % "Birth", 0th order reaction
        temp = find(A(1,:)==0);             % Randomly choose agent which becomes "alive" 
        A(2,temp(ceil(rand * size(temp,2)))) = 1;
    end
    
    tempA1 = find(A(1,:)==1);
    for i = 1:size(tempA1,2)
        if rand < P_d                       % "Death", 1st order reaction
            A(2,tempA1(i)) = 0;             % A agent dies
        end
    end
    
    Sa(t) = sum(A(2,:));
    A(1,:) = A(2,:);
    time(t) = time(t-1) + dt;
end

figure('Name','Birth-death process Time course','NumberTitle','off'); 
scatter(time,Sa,3,'xc');                                                hold on;
tmax = time(end);
[t, y_sol] = ode45(@bd_dif,0:tmax/500:tmax,Sa(1));
scatter(t,y_sol(:,1),3,'.b');                                           hold off;

axis([0 tmax 0 agents]);         
xlabel('t (sec)');                  ylabel('N_A(t)');    
legend('A stoc','A deter');

clear i t temp tempA1;
toc
