% This is the same as the Sigmoidal Response 2-step motif,
% except the steps occur serially, as would be expected in 
% a more realistic scenario.

% R <--> Rp <--> Rpp             Rates: k_f, k_r
%     ^       ^
%     |       |
%     S       S                  S: Signal (reactant in forward reactions)

% Assume number of S molecules does not change.
% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;             tic;
clc; 
%% Declare variables and functions
global k_f k_r S;

maxTime = 100;                          % Maximum Simulation time (sec)
dt = 1/100;                             % Constant (fixed) time step increment (sec)
t_steps = maxTime / dt;

agents = 100;          	disp(['Agents = ' num2str(agents)]);
k_f = 0.01;                            % basal forward rate (1st order)
k_r = 0.1;                             % basal reverse rate (1st order)
S = 40;                                 % Number of S molecules is constant

%% Initialize - Preallocate memory for variables; ** Initial conditions **
t = 1;                               % Time counter variable
Tr   = zeros(1,t_steps);             % Total sum of R agents in each time step	
Trp  = zeros(1,t_steps);             % Total sum of Rp agents in each time step
Trpp = zeros(1,t_steps);             % Total sum of Rpp agents in each time step

% ******** Initial conditions - Number of R, Rp, Rpp Agents ********
Tr(1) = floor(agents/4);      		Trp(1) = agents - Tr(1);	
Trpp(1) = 0;
% ******************************************************************

tempR = zeros(1,agents);                          
% Put "1" where agents are "alive", then randomize the array
for c=1:Tr(1),                      tempR(c)=1;             end
tempR = RandArray(tempR);           % Randomize array
Rv = [tempR ; tempR];               % Initialize vector storing state of R agents     
Rpv = ~Rv;                          % Initialize vector storing state of Rp agents
Rppv = zeros(2,agents);             % Initialize vector storing state of Rpp agents
% Notes on Rv, Rpv, Rppv:
% - Markov process, so only previous and current time steps needed --> 2 rows:  
% - R + Rp + Rpp = agents
clear c tempR;

%% ABK Simulation

P_f = 1 - exp(-k_f * S * dt);    % P_ber for forward reaction: 2nd order (but S is constant)         
% P_f = k_f * S * dt;              % P_dif for forward reaction: 2nd order
P_r = 1 - exp(-k_r * dt);        % P_ber for reverse reaction: 1st order            
% P_r = k_r * dt;                  % P_dif for reverse reaction: 1st order

% while abs(Trp(t)-Rp_ss) > 0.5 && t*dt <= maxTime
while t*dt <= maxTime
    
    for i = 1:agents
        if Rv(1,i) == 1                     % if R agent is alive
            if rand < P_f                  % **check probability condition**
                Rv(2,i) = 0;                 % R agent "dies"
                Rpv(2,i) = 1;                 % R is converted to Rp
            end
            
        elseif Rpv(1,i) == 1					% if Rp agent is "alive"
            r = rand;
            if r < P_r                  % **check probability condition**
                Rpv(2,i) = 0;                 % Rp agent "dies"
                Rv(2,i) = 1;                 % Rp is converted to R
            elseif r >= P_r && r < P_r + P_f
                Rpv(2,i) = 0;                 % Rp agent "dies"
                Rppv(2,i) = 1;                 % Rp is converted to Rpp
            end
 
        elseif Rppv(1,i) == 1					% if Rpp agent is "alive"
            if rand < P_r                  % **check probability condition**
                Rppv(2,i) = 0;                 % Rpp agent "dies"
                Rpv(2,i) = 1;                 % Rpp is converted to Rp
            end
        end
    end
    
    Tr(t+1) = sum(Rv(2,:));      Trp(t+1) = sum(Rpv(2,:));      Trpp(t+1) = sum(Rppv(2,:));
    
    Rv(1,:) = Rv(2,:);           Rpv(1,:) = Rpv(2,:);           Rppv(1,:) = Rppv(2,:); 
    
    t = t + 1;   
end

clear Rv Rpv Rppv;
finaltime = (t-1) * dt;

%% Solve differential equation (was waiting for variable 'finaltime')
[t_sol, y_sol] = ode45(@sigm_2step_real_dif,0:finaltime/100:finaltime+(finaltime/20),...
    [Tr(1) ; Trp(1) ; Trpp(1)]);
disp(['Diff eq terminal R  value   = ' num2str(y_sol(end,1))]);
disp(['Diff eq terminal Rp value   = ' num2str(y_sol(end,2))]);
disp(['Diff eq terminal Rpp value  = ' num2str(y_sol(end,3))]);

% Note: extended time domain of solution by finaltime/20 to get closer to steady-state values

%% Graph ABK (stochastic) and deterministic results
figure('Name','Time course','NumberTitle','off');
time = 0:dt:finaltime;
scatter(time,Tr,3,'oc');                                        hold on;
scatter(time,Trp,3,'om');
scatter(time,Trpp,3,'ok');

scatter(t_sol,y_sol(:,1),3,'b');                  
scatter(t_sol,y_sol(:,2),3,'r');  
scatter(t_sol,y_sol(:,3),3,'g');                               hold off;

axis([0 t_sol(end) 0 agents]);          xlabel('time');         ylabel('# Agents');         
legend('R stoc','Rp stoc','Rpp stoc','R deter','Rp deter','Rpp deter');

%% Finish
clear r w i j k temp;
toc
