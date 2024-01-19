%   Ep  <--->  E            Rates: k_f, k_r (Reverse rx activated by R)
%    |   ^                  Michaelis-Menten constants: Km_f, Km_r
%    |    \
%    \    |
%    -->  R   -->                      R: Response
%   k_b   ^   k_d
%         | k_s
%         S                            S: Signal

% Simulating mutual activation switch process for R: 
% k_b: R synthesis, 1st order process wrt Ep (but Ep is not consummed)
% k_b: ALTERNATIVE IMPLEMENTATION: R synthesis, 0th order process, UPregulated by Ep
% k_d: R degradation, 1st order process wrt R
% k_s: R synthesis, 1st order process wrt S
% k_s: ALTERNATIVE IMPLEMENTATION: R synthesis, 0th order process, UPregulated by S
% k_f: E synthesis, MM process
% k_r: Ep synthesis, MM process wrt R (but R is not consummed)
% Km_f, Km_r: MM constants for forward and reverse rxs

% Note: Assume number of S molecules/agents is NOT changing (eg, S is DNA molecules)

% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear;      clc;     tic;                           
rng(1);
reps = 100;                 % Number of times to repeat experiment
global agents k_b k_d k_s k_f k_r Km_f Km_r S;

totalTime = 1000;             % Simulation time (sec)
dt = 1/50;                   % Constant (fixed) time step increment (sec)
t_steps = totalTime / dt;
% time = zeros(1,t_steps);

agents = 100;
k_b = 0.02;                        % 0th order R synthesis rate, UPregulated by Ep
k_d = 0.075;                       % 1st order degradation rate (units: 1/sec)
k_s = 0.05;                        % 0th order R synthesis rate UPregulated by S
k_f = 1;                           % basal forward rate (1st order)
k_r = 0.05;                        % basal reverse rate (1st order)
Km_f = 5;                          % MICROSCOPIC Michaelis-Menten constant for forward rx
Km_r = 10;                         % MICROSCOPIC Michaelis-Menten constant for reverse rx

S = 15;                          % Assume number of S molecules/agents is NOT changing

% Initial population sizes
Ri = 10;        Epi = 50;      Ei = agents - Epi;					


[R_ss, Ep_ss] = RateCurve;
% figure('Name','Time Course','NumberTitle','off');     

for n = 1:reps
    
    Tr = zeros(1,t_steps);
    Tep = zeros(1,t_steps);             Te = zeros(1,t_steps);
    P_f = zeros(1,t_steps);             P_r = zeros(1,t_steps);
    P_b = zeros(1,t_steps);

    % ******** Initial conditions - Number of R, Ep, E Agents ********
    Tr(1) = Ri;             Tep(1) = Epi;           Te(1) = Ei;					
    % ****************************************************************

    tempR = zeros(1,agents);        tempEp = zeros(1,agents);             
    % Put "1" where agents are "alive", then randomize the array
    for c=1:Tr(1),                  tempR(c)=1;             end
    for d=1:Tep(1),                 tempEp(d)=1;            end
    tempR = RandArray(tempR);                   % Randomize R array
    tempEp = RandArray(tempEp);                 % Randomize Ep array
    % Markov process, so only previous and current time steps needed --> 2 rows:    
    R = [tempR ; tempR];                        % Initialize vector storing state of R agents
    Ep = [tempEp ; tempEp];                     % Initialize vector storing state of Ep agents
    E = ~ Ep;                                   % Initialize vector storing state of E agents
    % - Ep and E are complementary (E + Ep = agents)
    clear c d tempR tempEp;

    t = 1;

    P_d = k_d * dt;                 % Probability of R degradation process (1st order wrt R)
    P_s = k_s *S* dt;               % Probability of R synthesis process (0th order)

    while t <= t_steps % && Tep(t)<agents
        P_b(t) = k_b * Tep(t) * dt;      % Probability of R synthesis process (0th order)
        P_f(t) = k_f / (Km_f + Tep(t)) * dt;          % wrt each Ep molecule
        P_r(t) = k_r * Tr(t) / (Km_r + Te(t)) * dt;   % wrt each E molecule

        % Take care of 0th order processes first
        if rand < P_b(t)
            tempR = find(R(1,:)==0);    % Randomly choose R agent synthesis 
            R(2,tempR(ceil(rand * size(tempR,2)))) = 1;
        end

        if rand < P_s
            tempR = find(R(1,:)==0);    % Randomly choose R agent synthesis 
            R(2,tempR(ceil(rand * size(tempR,2)))) = 1;
        end
        % End of 0th order processes

        tempR = find(R(1,:)==1);
        for i = 1:size(tempR,2)
            if rand < P_d                   % Degradation of R, 1st order rx
                R(2,tempR(i)) = 0;          % R agent is degraded
            end
        end

        tempE = find(E(1,:)==1);
        for j=1:size(tempE,2)
            if rand < P_r(t)
                E(2,tempE(j)) = 0;          % Conversion of E to Ep
                tempEp = find(Ep(1,:)==0);  % Randomly choose Ep agent synthesis
                Ep(2,tempEp(ceil(rand * size(tempEp,2)))) = 1;
            end
        end

        tempEp = find(Ep(1,:)==1);
        for k=1:size(tempEp,2)
            if rand < P_f(t)
                Ep(2,tempEp(k)) = 0;        % Conversion of Ep to E
                tempE = find(E(1,:)==0);    % Randomly choose E agent synthesis
                E(2,tempE(ceil(rand * size(tempE,2)))) = 1;
            end
        end

        Tr(t+1) = sum(R(2,:));          
        Tep(t+1) = sum(Ep(2,:));            Te(t+1) = sum(E(2,:));
        R(1,:) = R(2,:);
        Ep(1,:) = Ep(2,:);                  E(1,:) = E(2,:);
        t = t + 1;
    end

    Rv(n,:) = Tr;                           Epv(n,:) = Tep;

    ABK_R_ss(n) = mean(Tr(600/dt:end));        % Avg R pop size from t=600-->totalTime
    
    finaltime = (t-1) * dt;

end                 % end 'for n'

%% Solve DE
[t_sol, y_sol] = ode45(@mAct_dif,0:totalTime/500:totalTime,[Ri; Epi; Ei]);
%% Plot (selected) time trajectories
time = 0:dt:finaltime;
scatter(t_sol,y_sol(:,1),3,'xc');                                       hold on;
scatter(t_sol,y_sol(:,2),3,'.m');
% scatter(t_sol,y_sol(:,3),3,'.g'); 
scatter(time,Rv(9,:),3,'.b');                                                
scatter(time,Epv(9,:),3,'or');
% scatter(time,Te,3,'og');
axis([0 finaltime 0 agents]);         xlabel('time');            ylabel('#');     hold off; 
% legend('R stoc','Ep stoc','E stoc','R deter','Ep deter','E deter','Location','Best');

%% Classify trajectories in terms of which equilibrium point they reach
sink1 = [];                sink2 = [];              % source = [];
for n=1:reps
    if ABK_R_ss(n) < 20
        sink1 = [sink1 n];
    elseif ABK_R_ss(n) > 30 % && Rv(n,end) < 600
        sink2 = [sink2 n];
    end
end
clear n;

%% Plot State Space - Modify accordingly
figure('Name','State Space','NumberTitle','off');               hold on;
scatter(y_sol(:,1),y_sol(:,2),'+k');                 % plot deterministic trajectory

% Modify the following if plotting a series of trajectories
% for n=1:size(sink1,2)
%     scatter(R(sink1(n),:),Ep(sink1(n),:),1,'.r');                    % plot stochastic trajectory
% end

% Modify the following for plotting specific trajectories
scatter(Rv(9,:),Epv(9,:),3,'.r');                    % plot stochastic trajectory
% scatter(Rv(39,:),Epv(39,:),3,'.g');                    % plot stochastic trajectory
% scatter(Rv(99,:),Epv(99,:),3,'.b');                    % plot stochastic trajectory

plot(R_ss,Ep_ss,'.c','MarkerSize',6);
xlabel('R');                    ylabel('E_p');                  hold off; 
legend('Deterministic','ABK stochastic');

%% Finish - Notes
clear h i j k n t Tr Tep Te R Ep E temp tempR tempEp tempE;
toc

% Notes:
% - Some trajectories don't follow the deterministic path and go to the other sink.
% - When R(1)=10, Ep(1)=50, S = 15, DE predicts converegence at fixed point at R=33.
% Approx. 45% of the time the trajectory went to the other sink, R=11.
