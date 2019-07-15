%    R  <--->  Rp            Rates: k_f, k_r (Forward rx activated by X)
%          ^    |            Michaelis-Menten constants: Km_f, Km_r
%         /     |
%         |     \
%    -->  X    -->                      Rp: Response
%   k_b   ^    k_d1/k_d2
%         | k_s
%         S                             S: Signal

% Simulating negative feedback process: 
% k_b: X synthesis, 0th order process
% k_s: X synthesis, 1st order process wrt S (but S is not consummed)
% k_s: ALTERNATIVE IMPLEMENTATION: X synthesis, 0th order process, UPregulated by S
% k_d1: X degradation, 1st order process wrt X
% k_d2: X degradation, 2nd order process wrt X, Rp (but Rp is not consummed)
% k_f: Rp synthesis, MM process wrt X (but X is not consummed)
% k_r: R synthesis, MM process 
% Km_f, Km_r: MM constants for forward and reverse rxs

% ** We simulate a system with delay in the negative feedback rx. **

% Note: Assume number of S molecules/agents is NOT changing (eg, S is DNA molecules)

% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;      clc;     tic;                           % rng(1);

global agents k_s k_b k_d1 k_d2 k_f k_r Km_f Km_r S;

totalTime = 5000;             % Simulation time (sec)
dt = 1/100;                   % Constant (fixed) time step increment (sec)
t_steps = totalTime / dt;
% time = zeros(1,t_steps);
lag = 50;                     % Explicit delay value for feedback reaction

agents = 100;
k_b = 0;                      % 0th order X synthesis rate
k_d1 = 0;                     % 1st order degradation rate (units: 1/sec)
k_d2 = 0.005;                 % 2nd order degradation rate wrt X, Rp (units: 1/sec)
k_s = 0.10;                   % 0th order R synthesis rate UPregulated by S
k_f = 0.01;                   % basal forward rate (1st order)
k_r = 1;                      % basal reverse rate (1st order)
Km_f = 10;                    % MICROSCOPIC Michaelis-Menten constant for forward rx
Km_r = 10;                    % MICROSCOPIC Michaelis-Menten constant for reverse rx

S = 30;                       % Assume number of S molecules/agents is NOT changing

Xi = 0;
Rpi = 0;
%% Initialize, set Initial conditions.
Tx = zeros(1,t_steps);
Tr = zeros(1,t_steps);              Trp = zeros(1,t_steps);
P_f = zeros(1,t_steps);             P_r = zeros(1,t_steps);
P_d2 = zeros(1,t_steps);

% ******** Initial conditions - Number of R, Ep, E Agents ********
Tx(1) = Xi;         Trp(1) = Rpi;         Tr(1) = agents - Trp(1);					
% ****************************************************************

tempX = zeros(1,2*agents);        tempR = zeros(1,agents);             
% Put "1" where agents are "alive", then randomize the array
for c=1:Tx(1),                  tempX(c)=1;             end
for d=1:Tr(1),                  tempR(d)=1;            end
tempX = RandArray(tempX);                   % Randomize R array
tempR = RandArray(tempR);                   % Randomize R array
% Markov process, so only previous and current time steps needed --> 2 rows:    
Xv = [tempX ; tempX];                        % Initialize Vector storing state of X agents
Rv = [tempR ; tempR];                        % Initialize Vector storing state of R agents
Rpv = ~ Rv;                                  % Initialize Vector storing state of Rp agents
% - R and Rp are complementary (R + Rp = agents)
clear c d tempX tempR;
%% ABK simulation
t = 1;

P_b = k_b * dt;                  % Probability of X synthesis process (0th order)
P_s = k_s * S * dt;              % Probability of X synthesis process (0th order), UPreg'd by S
P_d1 = 1-exp(-k_d1 * dt);        % Probability of X degradation process (1st order) wrt X

while t <= t_steps % && Tep(t)<agents
    
    % ****** Introducing Delay in the Probability Expression of Feedback Rx ******     
    if t <= lag/dt
        P_d2(t) = k_d2 * Trp(t) * dt;                      % wrt each X molecule
    else
        P_d2(t) = k_d2 * Trp(t-lag/dt) * dt;           % introduce delay of 50 sec
    end
    % ****************************************************************************
    
    P_f(t) = k_f * Tx(t) / (Km_f + Tr(t)) * dt;            % wrt each R molecule
    P_r(t) = k_r / (Km_r + Trp(t)) * dt;                   % wrt each Rp molecule
    
    % Take care of 0th order processes first
    if rand < P_b
        tempX = find(Xv(1,:)==0);                % Randomly choose X agent synthesis 
        Xv(2,tempX(ceil(rand * size(tempX,2)))) = 1;
    end
    
    if rand < P_s
        tempX = find(Xv(1,:)==0);                % Randomly choose X agent synthesis 
        Xv(2,tempX(ceil(rand * size(tempX,2)))) = 1;
    end
    % End of 0th order processes
    
% The following implementation also works (treating it as a 1st order rx wrt S)     
%     for h=1:S                               % Reaction for each S molecule/agent
%         if rand < P_s
%             tempX = find(Xv(1,:)==0);       % Randomly choose X agent which becomes "alive" 
%             Xv(2,tempR(ceil(rand * size(tempX,2)))) = 1;
%         end
%     end
    
    P_d = P_d1 + P_d2(t);               % Total probability of X degradation
    tempX = find(Xv(1,:)==1);
    for i = 1:size(tempX,2)
        if rand < P_d                   % Degradation of X, 2 contributing processes
            Xv(2,tempX(i)) = 0;         % X agent is degraded
        end
    end

    tempR = find(Rv(1,:)==1);
    for j=1:size(tempR,2)
        if rand < P_f(t)
            Rv(2,tempR(j)) = 0;         % Conversion of R ...
            Rpv(2,tempR(j)) = 1;        % to Rp
        end
    end
    
    tempRp = find(Rpv(1,:)==1);
    for k=1:size(tempRp,2)
        if rand < P_r(t)
            Rpv(2,tempRp(k)) = 0;       % Conversion of Rp ...
            Rv(2,tempRp(k)) = 1;        % to R
        end
    end
    
    Tx(t+1) = sum(Xv(2,:));
    Tr(t+1) = sum(Rv(2,:));             Trp(t+1) = sum(Rpv(2,:));      
 
    Xv(1,:) = Xv(2,:);
    Rv(1,:) = Rv(2,:);                  Rpv(1,:) = Rpv(2,:);
    t = t + 1;
end

% Remove unnecessary terminal 0's from arrays
if t < t_steps
    Tx = Tx(1:t);               Tr = Tr(1:t);           Trp = Trp(1:t);                                           
    P_f = P_f(1:t);             P_r = P_r(1:t);         P_d2 = P_d2(1:t);                 
end

%% Solve Delay DE 
global DEhistory;
% Solve ODE for the initial time period 0:lag
DEhistory = ode23(@negFb_dif,0:dt:lag,[Xi; Rpi]);   % DEhistory is a structure

cd DDE_files;
DDEsol = dde23(@negFb_delaydif,lag,@negFb2_dde_hist,[lag , totalTime]);

% plot(DDEsol.x,DDEsol.y);

t_sol = [lag:totalTime]';
y_sol = deval(DDEsol,t_sol)';                           cd ..;
%% Plot time course
time = 0:dt:totalTime;
figure('Name','Time Course','NumberTitle','off');               hold on;
plot(time,Tx,'b');                                
plot(time,Trp,'r');
plot(t_sol,y_sol(:,1),'c');                                       
plot(t_sol,y_sol(:,2),'b');
axis([0 totalTime 0 1.5*agents]);                               hold off;
xlabel('time');                 ylabel('#');      
legend('X stoc','Rp stoc','X deter','Rp deter','Location','Best');

%% Finish
clear array X Rp dX dRp X_nc Rp_nc X_nc_sym Rp_nc_sym;
clear h i j k n r t temp tempR tempRp tempX Xv Rv Rpv;
toc
