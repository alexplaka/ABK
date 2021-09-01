%   Ep <--->  E            Rates: k_f, k_r (Reverse rx activated by R)
%    |  ^                  Michaelis-Menten constants: Km_f, Km_r
%    |   \   --> X -->     k_bx and k_dx: synthesis and degradation of X
%    \   \  /    \
%    -->  R  -------->                      R: Response
%   k_b   ^      k_d
%         | k_s
%         S                            S: Signal

% Simulating mutual activation switch process for R: 
% k_b: R synthesis, 1st order process wrt Ep (but Ep is not consummed)
% k_d1: R degradation, 1st order process wrt R
% k_d2: R degradation, 2nd order process wrt R, X
% k_bx: X synthesis, 0th order, UPregulated by R
% k_dx: X degradation, 1st order wrt X
% k_s: R synthesis, 1st order process wrt S
% k_f: E synthesis, MM process
% k_r: Ep synthesis, MM process wrt R (but R is not consummed)
% Km_f, Km_r: MM constants for forward and reverse rxs

% Note: Assume number of S molecules/agents is NOT changing (eg, S is DNA molecules)

% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;      clc;     tic;
rng(0);
global agents k_b k_d1 k_d2 k_s k_bx k_dx k_f k_r Km_f Km_r S;

totalTime = 2000;             % Simulation time (sec)
dt = 1/50;                    % Constant (fixed) time step increment (sec)
t_steps = totalTime / dt;
% time = zeros(1,t_steps);

agents = 100;
k_b = 0.01;                        % 0th order R synthesis rate, UPregulated by Ep
k_d1 = 0.0;                        % 1st order degradation rate (units: 1/sec)
k_d2 = 0.0075;                     % MICROscopic 2nd order degradation rate (units: 1/sec)
k_s = 0.05;                        % 0th order R synthesis rate UPregulated by S
k_bx = 0.001;                      % 0th order synthesis of X, UPregulated by R
k_dx = 0.01;                       % 1st order degradation of X
k_f = 1;                           % basal forward rate (1st order)
k_r = 0.05;                        % basal reverse rate (1st order)
Km_f = 5;                          % MICROSCOPIC Michaelis-Menten constant for forward rx
Km_r = 10;                         % MICROSCOPIC Michaelis-Menten constant for reverse rx
S = 15;                            % Assume number of S molecules/agents is NOT changing

% Initial population sizes
Ri = 0;                 Epi = 0;                Xi = 0;

% Initialize - Preallocate memory
Tr = zeros(1,t_steps);              Tx = zeros(1,t_steps);
Tep = zeros(1,t_steps);             Te = zeros(1,t_steps);

P_f = zeros(1,t_steps);             P_r = zeros(1,t_steps);
P_b = zeros(1,t_steps);             P_d2 = zeros(1,t_steps);
P_bx = zeros(1,t_steps);

% ******** Initial conditions - Number of R, Ep, E Agents ********
Tr(1) = Ri;                       Tx(1) = Xi;    
Tep(1) = Epi;                     Te(1) = agents - Epi;					
% ****************************************************************

tempR = zeros(1,2*agents);        tempEp = zeros(1,agents); 
tempX = zeros(1,agents);
% Put "1" where agents are "alive", then randomize the array
for c=1:Tr(1),                  tempR(c)=1;             end
for d=1:Tep(1),                 tempEp(d)=1;            end
for e=1:Tx(1),                  tempX(e)=1;             end
tempR = RandArray(tempR);                   % Randomize R array
tempEp = RandArray(tempEp);                 % Randomize Ep array
tempX = RandArray(tempX);                   % Randomize X array
% Markov process, so only previous and current time steps needed --> 2 rows:    
R = [tempR ; tempR];                        % Initialize vector storing state of R agents
Ep = [tempEp ; tempEp];                     % Initialize vector storing state of Ep agents
E = ~ Ep;                                   % Initialize vector storing state of E agents
% - Ep and E are complementary (E + Ep = agents)
X = [tempX ; tempX];                        % Initialize vector storing state of X agents
clear c d e tempR tempEp tempX;

% [R_ss, Ep_ss] = RateCurve;           disp(['All R_ss values = ' num2str(R_ss')]);
% bifurc(Tr(1),Tep(1));
%% ABK simulation
t = 1;

P_s = k_s * S * dt;             % Probability of R synthesis process (0th order)
P_d1 = k_d1 * dt;               % Probability of R degradation (1st order)
P_dx = k_dx * dt;               % Probability of X degradation process (1st order wrt X)

while t <= t_steps % && Tep(t)<agents
    
    P_b(t) = k_b * Tep(t) * dt;      % Probability of R synthesis process (0th order)
    P_bx(t)= k_bx * Tr(t) * dt;     % Probability of X synthesis process (0th order)
    P_f(t) = k_f / (Km_f + Tep(t)) * dt;          % wrt each Ep molecule
    P_r(t) = k_r * Tr(t) / (Km_r + Te(t)) * dt;   % wrt each E molecule
    P_d2(t) = k_d2 * Tx(t) * dt;                   % wrt each R molecule
    
    % Take care of 0th order processes first
    if rand < P_b(t)
        tempR = find(R(1,:)==0);    % Randomly choose R agent synthesis 
        R(2,tempR(ceil(rand * size(tempR,2)))) = 1;
    end

    if rand < P_bx(t)
        tempX = find(X(1,:)==0);    % Randomly choose X agent synthesis 
        X(2,tempX(ceil(rand * size(tempX,2)))) = 1;
    end
    
    if rand < P_s
        tempR = find(R(1,:)==0);    % Randomly choose R agent synthesis 
        R(2,tempR(ceil(rand * size(tempR,2)))) = 1;
    end
    % End of 0th order processes
    
% The following implementation also works (treating it as a 1st order rx wrt S)     
%     for h=1:S                               % Reaction for each S molecule/agent
%         if rand < P_s
%             tempR = find(R(1,:)==0);         % Randomly choose R agent which becomes "alive" 
%             R(2,tempR(ceil(rand * size(tempR,2)))) = 1;
%         end
%     end
    
    P_d = P_d1 + P_d2(t);               % Total probability of R degradation
    tempR = find(R(1,:)==1);
    for i = 1:size(tempR,2)
        if rand < P_d                   % Degradation of R, 1st + 2nd order rx
            R(2,tempR(i)) = 0;          % R agent is degraded
        end
    end

    tempE = find(E(1,:)==1);
    for j=1:size(tempE,2)
        if rand < P_r(t)
            E(2,tempE(j)) = 0;          % Conversion of E to Ep
            Ep(2,tempE(j)) = 1;
        end
    end
    
    tempEp = find(Ep(1,:)==1);
    for k=1:size(tempEp,2)
        if rand < P_f(t)
            Ep(2,tempEp(k)) = 0;        % Conversion of Ep to E
            E(2,tempEp(k)) = 1;
        end
    end
    
    tempX = find(X(1,:)==1);
    for m = 1:size(tempX,2)
        if rand < P_dx                  % Degradation of X, 1st order rx
            X(2,tempX(m)) = 0;          % X agent is degraded
        end
    end  
    
    Tr(t+1) = sum(R(2,:));              Tx(t+1) = sum(X(2,:));          
    Tep(t+1) = sum(Ep(2,:));            Te(t+1) = sum(E(2,:));
    
    R(1,:) = R(2,:);                    X(1,:) = X(2,:);
    Ep(1,:) = Ep(2,:);                  E(1,:) = E(2,:);
    
    t = t + 1;
end

% Remove unnecessary terminal 0's from arrays
if t < t_steps
    Tr = Tr(1:t);                               Tx = Tx(1:t);                             
    Tep = Tep(1:t);                             Te = Te(1:t);             
    P_f = P_f(1:t);                             P_r = P_r(1:t);
    P_d2 = P_d2(1:t);      P_b = P_b(1:t);      P_bx = P_bx(1:t);
end

clear R Ep E X;
%% Solve DE
if exist('t','var')==0,         finaltime = totalTime;      
else                            finaltime = (t-1) * dt;         end
[t_sol, y_sol] = ode45(@actinh2_dif,0:finaltime/500:finaltime,[Ri ; Epi ; Xi]);
%% Plot time course
time = 0:dt:finaltime;
figure('Name','Time Course','NumberTitle','off');                           hold on;

plot(time(1:20:end),Tr(1:20:end),'r');                                
plot(time(1:20:end),Tep(1:20:end),'g');
plot(time(1:20:end),Tx(1:20:end),'b');
plot(t_sol,y_sol(:,1),'--','Color',[0.75 0 1],'LineWidth',2);                                       
plot(t_sol,y_sol(:,2),'--','Color',[0 0.5 0],'LineWidth',2);
plot(t_sol,y_sol(:,3),'--','Color',[0.15 1 0.75],'LineWidth',2);                                     

axis([0 finaltime 0 agents]);
% axis([0 finaltime 0 max(Tr)]);                               
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');                   hold off;
xlabel('t (sec)');                 ylabel('N(t)');      
leg = legend('ABK N_R(t)','ABK N_{Ep}(t)','ABK N_X(t)','DE N_R(t)','DE N_{Ep}(t)','DE N_X(t)');
set(leg,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95],'Location','Best');                                              

%% Finish
clear array R Ep dR dEp R_nc Ep_nc R_nc_sym Ep_nc_sym;
clear h i j k n r t temp tempR tempEp tempE leg;
toc
