% This is the same as the Hyperbolic Response system,
% except the forward rx depends on S^2, as would be 
% the case if two molecules of S were needed to 
% accomplish the conversion R --> Rp.

% R <--> Rp             Rates: k_f, k_r
%     ^
%     |
%     2S                S: Signal (reactant in forward reaction)

% Assume number of S molecules does NOT change.

% Considering HETEROGENEITY in Rp agents' kinetic rate constant.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear;             tic;                     clc; 
rng(0);
%% Declare variables and functions
global k_f k_r_mean S;

maxTime = 100;                          % Maximum Simulation time (sec)
dt = 1/100;                             % Constant (fixed) time step increment (sec)
t_steps = maxTime / dt;

agents = 100;          	disp(['Agents = ' num2str(agents)]);

k_f = 0.0001;                          % basal forward rate (1st order)

k_r_mean = 0.1;                        % mean basal reverse rate (1st order)
k_r_sdev = 0.025;                      % sdev of basal reverse rate (1st order)
k_r = [k_r_mean + k_r_sdev * randn(1,agents)]; % PHM for k_r values

S = 7;                                % Number of S molecules is constant
 
%% Simple formula for Rp_ss
Rp_ss = agents * S^2 / ( k_r_mean/k_f + S^2 );

%% Initialize - Preallocate memory for variables; ** Initial conditions **
t = 1;                              % Time counter variable
Tr = zeros(1,t_steps);              % Total sum of R agents in each time step	
Trp = zeros(1,t_steps);             % Total sum of Rp agents in each time step

% ******** Initial conditions - Number of R, Rp Agents ********
Tr(1) = floor(agents/4);      		Trp(1) = agents - Tr(1);					
% ************************************************************

tempR = zeros(1,agents);                          
% Put "1" where agenst are "alive", then randomize the array
for c=1:Tr(1),                      tempR(c)=1;             end
tempR = RandArray(tempR);           % Randomize array
Rv = [tempR ; tempR];               % Initialize vector storing state of R agents     
Rpv = ~Rv;                          % Initialize vector storing state of Rp agents 
% Notes on Rv, Rpv:
% - Markov process, so only previous and current time steps needed --> 2 rows:          
% - Rv and Rpv are complementary (R+Rp=agents)

avg_k_r = zeros(1,t_steps+1);       % For storing instantaneous <k_r(t)>
tempRp_log = Rpv(2,:)==1;           % Logical array of existing E agents at end of time step
avg_k_r(1) = tempRp_log * k_r' / sum(tempRp_log); 

clear c tempR*;

%% ABK Simulation

P_f = 1 - exp(-k_f * S^2 * dt);    % P_ber for forward reaction: 1st order            
%     P_f = k_f * S^2 * dt;              % P_dif for forward reaction: 1st order
P_r = 1 - exp(-k_r * dt);        % P_ber for reverse reaction: 1st order  (using k_r PHM)        
%     P_r = k_r * dt;                  % P_dif for reverse reactio: 1st order (using k_r PHM)   

% while abs(Trp(t)-Rp_ss) > 0.5 && t*dt <= maxTime
while t*dt <= maxTime
    
    for i = 1:agents
        if Rv(1,i) == 1                     % if R agent is alive
            if rand < P_f                  % **check probability condition**
                Rv(2,i) = 0;                 % R agent "dies"
                Rpv(2,i) = 1;                 % R is converted to Rp
            end
        elseif Rpv(1,i) == 1					% if Rp agent is "alive"
            if rand < P_r(i)                  % **check probability condition**
                Rpv(2,i) = 0;                 % Rp agent "dies"
                Rv(2,i) = 1;                 % Rp is converted to R
            end
        end
    end
    
    Tr(t+1) = sum(Rv(2,:));                  Trp(t+1) = sum(Rpv(2,:));
    Rv(1,:) = Rv(2,:);                       Rpv(1,:) = Rpv(2,:); 
    
    tempRp_log = Rpv(2,:)==1;     % Logical array of existing Rp agents at end of time step
    avg_k_r(t+1) = tempRp_log * k_r' / sum(tempRp_log);    % Find <k_r> at end of time step

    t = t + 1;   
    
end

% Remove unnecessary terminal 0's from arrays
Tr = Tr(Tr~=0);                             w = size(Tr,2);                 
Trp = Trp(1:w);                                      
disp(['ABK sim terminal R value   = ' num2str(Tr(end))]);  
disp(['ABK sim terminal Rp value  = ' num2str(Trp(end))]);        
clear Rv Rpv;

finaltime = (t-1) * dt;
time = 0:dt:finaltime;
%% Solve differential equation (was waiting for variable 'finaltime')
[t_sol, y_sol] = ode45(@sigm_2step_dif,0:finaltime/100:finaltime+(finaltime/20),[Tr(1) ; Trp(1)]);
disp(['Diff eq terminal R  value   = ' num2str(y_sol(end,1))]);
disp(['Diff eq terminal Rp value   = ' num2str(y_sol(end,2))]);
% Note: extended time domain of solution by finaltime/20 to get closer to steady-state values

%% Graph ABK (stochastic) and deterministic results
figure('Name','Time course','NumberTitle','off');

scatter(time,Tr,3,'oc');                                        hold on;
scatter(time,Trp,3,'om');
scatter(t_sol,y_sol(:,1),3,'.b');                  
scatter(t_sol,y_sol(:,2),3,'.r');                               hold off;
axis([0 t_sol(end) 0 agents]);          xlabel('time');         ylabel('# Agents');         
legend('R stoc','Rp stoc','R deter','Rp deter');

%%
figure
plot(time,avg_k_r);
xlabel('time');         ylabel('<k_r>');         

%% Finish
clear r w i j k temp;
toc
