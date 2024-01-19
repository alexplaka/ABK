% Hyperbolic Response - Adaptation system 

%             |                        Rate: k_sp
%             \
% --> R <---> Rp -->                   Rates: k_sr, k_f, k_r, k_dp
%        /      /
%        |     /                       S: Signal (reactant in R forward reaction)
%        S    /             
%        |    |
%        \    |                        S: Signal (upregulates synthesis of X)
%        -->  X  -->               
%        k_sx    k_dx                  Synthesis and degradation rates of X

% k_sr:  0th order
% k_f:   2nd order wrt S, R
% k_r:   1st order wrt Rp
% k_sp:  0th order
% k_dp:  2nd order wrt Rp, X
% k_sx:  0th order (upregulated by S)
% k_dx:  1st order wrt X

% Assume number of S molecules does NOT change.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear;             tic;
clc; 
rng(0);
%% Declare variables and functions
global k_sr k_f k_r k_sp k_dp k_sx k_dx S;

N_avo = 6.02e23;        % Avogadro's number
V = 1e-21;              % Volume in L

maxTime = 300;                        % Maximum Simulation time (sec)
dt = 1/50;                              % Constant (fixed) time step increment (sec)
t_steps = maxTime / dt;

agents = 500;          	disp(['Agents = ' num2str(agents)]);

ko_sr = 0.05;                           % basal rate of R synthesis (0th order)
k_sr = ko_sr * (N_avo * V);
ko_f = 10;                               % basal forward rate (2nd order)
k_f = ko_f / (N_avo * V);
k_r = 0.01;                             % basal reverse rate (1st order)
ko_sp = 1e-4;                           % basal rate of Rp synthesis (0th order)
k_sp = ko_sp * (N_avo * V);
ko_dp = 2;                              % basal rate of Rp degradation (2nd order)
k_dp = ko_dp / (N_avo * V);
ko_sx = 1e-4;                           % basal rate of X synthesis (0th order)
k_sx = ko_sx * (N_avo * V);         
k_dx = 0.01;                           % basal rate of X degradation (1st order)

S = 30;                                 % Number of S molecules is constant

%% *** Set up system of differential equations using symbolic variables ***
% *** (requires Symbolic Math Toolbox) ***
syms R Rp X positive;
dR_sym  = + k_sr - k_f * R * S + k_r * Rp;
dRp_sym = + k_sp + k_f * R * S - k_r * Rp - k_dp * Rp * X;
dX_sym  = + k_sx * S - k_dx * X;

% Nullclines (Equilibrium conditions)
R_eq_sym = solve(dR_sym == 0,R);
Rp_eq_sym = solve(dRp_sym == 0,Rp);
X_eq_sym = solve(dX_sym == 0,X);

%% Symbolic calculations
% Calculate theoretical steady state values for R, Rp, X:
X_ss = double(X_eq_sym);
% Rp_eq_sub = subs(dRp_sym,[R X],[R_eq_sym X_ss]);
% Rp_ss = double(solve(Rp_eq_sub==0,Rp));
% R_ss = double(subs(R_eq_sym,Rp,Rp_ss));
% 
% disp(['Theoretical  R_ss = ' num2str(R_ss)]);
% disp(['Theoretical Rp_ss = ' num2str(Rp_ss)]);
disp(['Theoretical  X_ss = ' num2str(X_ss)]);

%% Initialize - Preallocate memory for variables; ** Initial conditions **
t = 1;                              % Time counter variable
time = zeros(1,t_steps);
P_dp = zeros(1,t_steps);            % For storing probability value at each time step

Tr = zeros(1,t_steps);              % Total sum of R agents in each time step	
Trp = zeros(1,t_steps);             % Total sum of Rp agents in each time step
Tx = zeros(1,t_steps);              % Total sum of X agents in each time step	

% ******** Initial conditions - Number of R, Rp Agents *****************************
Tr(1) = 0;                		Trp(1) = 0;               Tx(1) = 0;				
% **********************************************************************************
Rv = zeros(2,agents);           Rpv = zeros(2,agents);      Xv = zeros(2,agents);

% % If initial condition is not all zeros, then run below code instead.
% tempR = zeros(1,agents);                          
% % Put "1" where agenst are "alive", then randomize the array
% for c=1:Tr(1),                      tempR(c)=1;                 end
% tempR = RandArray(tempR);           % Randomize array
% Rv = [tempR ; tempR];               % Initialize vector storing state of R agents     
% Rpv = ~Rv;                          % Initialize vector storing state of Rp agents 
% % Notes on Rv, Rpv:
% % - Markov process, so only previous and current time steps needed --> 2 rows:          
% % - Rv and Rpv are complementary (R+Rp=agents) ONLY AS AN INITIAL CONDITION.
% 
% tempX = zeros(1,agents);             
% % Put "1" where agents are "alive", then randomize the array
% for d=1:Tx(1),                      tempX(d)=1;                 end
% tempX = RandArray(tempX);                   % Randomize array
% % Markov process, so only previous and current time steps needed --> 2 rows:    
% Xv = [tempX ; tempX];                       % Initialize vector storing state of X agents
% clear c d tempR tempX;

%% ABM Simulation
P_sr = k_sr * dt;
P_f =  k_f * S * dt;
P_r =  k_r * dt;
P_sp = k_sp * dt;
P_sx = k_sx * S * dt;
P_dx = k_dx * dt;

while t*dt <= maxTime && Tx(t)<agents
    P_dp(t) = k_dp * Tx(t) * dt;
	
    % Take care of 0th order processes (synthesis of R, X) first
    if rand < P_sr
        tempR = find(Rv(1,:)==0);                  % Randomly choose R agent which becomes "alive" 
        Rv(2,tempR(ceil(rand * size(tempR,2)))) = 1;
    end
    
    if rand < P_sp
        tempRp = find(Rpv(1,:)==0);                  % Randomly choose Rp agent which becomes "alive" 
        Rpv(2,tempRp(ceil(rand * size(tempRp,2)))) = 1;
    end
   
    if rand < P_sx
        tempX = find(Xv(1,:)==0);                  % Randomly choose X agent which becomes "alive" 
        Xv(2,tempX(ceil(rand * size(tempX,2)))) = 1;
    end
    % ** End of 0th order rxs **
    
    tempR = find(Rv(1,:)==1);
    for h=1:size(tempR,2)                          % Reaction for each R molecule/agent
        if rand < P_f
            Rv(2,tempR(h)) = 0;                    % Conversion of R to Rp
            tempRp = find(Rpv(1,:)==0);            % Randomly choose Rp agent which becomes "alive" 
            Rpv(2,tempRp(ceil(rand * size(tempRp,2)))) = 1;
        end
    end

    tempRp = find(Rpv(1,:)==1);
    for j=1:size(tempRp,2)                         % Reaction for each Rp molecule/agent
        r =  rand;
        if r < P_r
            Rpv(2,tempRp(j)) = 0;                  % Conversion of Rp to R
            tempR = find(Rv(1,:)==0);              % Randomly choose R agent which becomes "alive"
            Rv(2,tempR(ceil(rand * size(tempR,2)))) = 1;
        elseif r >= P_r && r < P_dp(t)
            Rpv(2,tempRp(j)) = 0;                  % Degradation of Rp
        end
    end
    
    tempX = find(Xv(1,:)==1);
    for q=1:size(tempX,2)                          % Reaction for each X molecule/agent
        if rand < P_dx
            Xv(2,tempX(q)) = 0;                    % Degradation of X
        end
    end
    
    time(t+1) = time(t) + dt;
    Tr(t+1) = sum(Rv(2,:));           Trp(t+1) = sum(Rpv(2,:));         Tx(t+1) = sum(Xv(2,:));
    Rv(1,:) = Rv(2,:);                Rpv(1,:) = Rpv(2,:);              Xv(1,:) = Xv(2,:);
    t = t + 1;   
end

% Remove unnecessary terminal 0's from arrays
if t < t_steps
    Tr = Tr(1:t);                                         
    Trp = Trp(1:t);                             Tx = Tx(1:t);             
    P_dp = P_dp(1:t);                           time = time(1:t);
end

disp(['ABK sim terminal R value   = ' num2str(Tr(end))]);  
disp(['ABK sim terminal Rp value  = ' num2str(Trp(end))]); 
disp(['ABK sim terminal X value   = ' num2str(Tx(end))]);                   clear Rv Rpv Xv;

finaltime = (t-1) * dt;
% time = 0:dt:finaltime;

%% Solve differential equation (was waiting for variable 'finaltime')
[t_sol, y_sol] = ode45(@hadapt_dif,0:finaltime/100:finaltime,[Tr(1) ; Trp(1); Tx(1)]);
disp(['Diff eq terminal R value   = ' num2str(y_sol(end,1))]);
disp(['Diff eq terminal Rp value  = ' num2str(y_sol(end,2))]);
disp(['Diff eq terminal X value   = ' num2str(y_sol(end,3))]);

%% Graph ABK (stochastic) and deterministic results
figure('Name','Time course','NumberTitle','off');
time = 0:dt:finaltime;
scatter(time,Tr,3,'oc');                                        hold on;
scatter(time,Trp,3,'ok');
scatter(time,Tx,3,'og');  
scatter(t_sol,y_sol(:,1),3,'xc');                  
scatter(t_sol,y_sol(:,2),3,'xk');
scatter(t_sol,y_sol(:,3),3,'xg');                               hold off;
axis([0 t_sol(end) 0 agents]);          xlabel('time');         ylabel('# Agents');         
legend('R stoc','Rp stoc','X stoc','R deter','Rp deter','X deter');

%% Finish
clear r w h j q t tempR tempRp tempX;
toc
