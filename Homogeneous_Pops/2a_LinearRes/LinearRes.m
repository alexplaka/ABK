%  k_b   k_d
%  --> R -->                        R: Response
%      ^
%      | k_s
%      S                            S: Signal

% Simulating birth-death process for R, and its synthesis from S: 
% k_b: birth, 0th order process; k_d: death, 1st order process
% k_s: synthesis from S, 1st order process
% Note: Assume number of S molecules/agents is NOT changing (eg, S is DNA molecules)
% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear;      clc;     tic;
rng(0);
agents = 1200;

N_avo = 6.02e23;                % Avogadro's number
V = 1e-21;                      % Volume in L

totalTime = 700;                % Simulation time (sec)
dt = 1/50;                      % Fixed time step increment (sec)

if exist('dt','var') == 1         
    t_steps = totalTime / dt;
else
    t_steps = 1200;             % empirically chosen value (for lambda = 0.8)
end
time = zeros(1,t_steps);

global k_b k_d k_s S;
k_bm = 3.322e-5;                % Molar "birth" rate; 0th order (units: M/sec)
k_b = k_bm * N_avo * V;         % Microscopic 0th order birth rate
k_d = 0.02;                     % "Death" rate; 1st order (units: 1/sec)
k_s = 1;                        % Synthesis rate; 1st order (units: 1/sec)

S = 20;                         % Do a simulation for this value of S
R_ss = (k_b + k_s * S) / k_d;   % Deterministic Steady-state value for R
    
Tr = zeros(1,t_steps);              Tr(1) = 0;        % ** Initial condition **

tempR = zeros(1,agents);             
% Put "1" where agents are "alive", then randomize the array
for c=1:Tr(1),                  tempR(c)=1;             end
tempR = RandArray(tempR);                   % Randomize array
% Markov process, so only previous and current time steps needed --> 2 rows:    
R = [tempR ; tempR];                        % Initialize vector storing state of A agents
clear c tempR;

P_b = k_b * dt;               % P_dif of "birth" process (0th order)

P_d = 1 -  exp(-k_d * dt);       % P_ber of "death" process (1st order wrt R)
% P_d = k_d * dt;               % P_dif of "death" process (1st order wrt R)

P_s0 = k_s * S * dt;          % P_dif of synthesis process (0th order wrt S)

P_s1 = 1 - exp(-k_s * dt);       % P_ber of synthesis process (1st order wrt S)
% P_s1 = k_s * dt;              % P_dif of synthesis process (1st order wrt S)

for t = 2:t_steps

    if rand < P_b                           % "Birth", 0th order reaction
        temp = find(R(1,:)==0);             % Randomly choose R agent which becomes "alive" 
        R(2,temp(ceil(rand * size(temp,2)))) = 1;
    end

%   --------------------------------------------------------------
%   Can treat production of R through S either as 0th or 1st order
    if rand < P_s0                          % "Birth", 0th order reaction promoted by S
        temp = find(R(1,:)==0);             % Randomly choose R agent which becomes "alive" 
        R(2,temp(ceil(rand * size(temp,2)))) = 1;
    end    

%     for h=1:S                               % Reaction for each S molecule/agent (1st order)
%         if rand < P_s1
%             temp = find(R(1,:)==0);         % Randomly choose R agent which becomes "alive" 
%             R(2,temp(ceil(rand * size(temp,2)))) = 1;
%         end
%     end
%   --------------------------------------------------------------

    for i = 1:size(R,2)
        if R(1,i) == 1                      % if R agent still exists          
            if rand < P_d                   % "Death", 1st order reaction
                R(2,i) = 0;                 % R agent is degraded
            end
        end
    end

    Tr(t) = sum(R(2,:));
    R(1,:) = R(2,:);
    time(t) = time(t-1) + dt;
end

%   Find average R of last 300sec (steady-state has been reached)
ABK_r = mean(Tr(end-300*50:end));

figure0 = figure('Name','Time Course','NumberTitle','off');
set(figure0,'Position',[1 1 500 406]); 
scatter(time,Tr,3,'xc');                                                hold on;
tmax = time(end);
[t, y_sol] = ode45(@linearRes_dif,0:tmax/500:tmax,Tr(1));
scatter(t,y_sol(:,1),3,'.b');                                           hold off;
axis([0 tmax 0 agents]);         
xlabel('t (sec)');            ylabel('N_R(t)'); 
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
leg0 = legend('R(t) stoc','R(t) deter');
set(leg0,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95]);

RateCurve(agents);              

%% Finish
clear h i t temp R;
toc
