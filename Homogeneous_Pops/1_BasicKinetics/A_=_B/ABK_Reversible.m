% A <--> B  
% Simulating 1st order kinetics for reversible chemical reaction.
% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;            tic;
clc; 
rng(0);
%% Declare variables and functions
global ko_f ko_r;

maxTime = 50;                        % Maximum Simulation time (sec)
dt = 1/100;                             % Constant (fixed) time step increment (sec)
t_steps = maxTime / dt;
agents = 100;          	disp(['Agents = ' num2str(agents)]);
ko_f = 0.02;                               % basal forward rate
ko_r = 0.01;                               % basal reverse rate

%% Initialize - Preallocate memory for variables; ** Initial conditions **
t = 1;                              % Time counter variable
P = zeros(1,t_steps);               % For storing probability value at each time step
Sa = zeros(1,t_steps);              % Sum of A agents in each time step	
Sb = zeros(1,t_steps);              % Sum of B agents in each time step

% ******** Initial conditions - Number of A, B Agents ********
Sa(1) = floor(agents/2);            Sb(1) = agents - Sa(1);					
% ************************************************************

tempA = zeros(1,agents);            % tempB = zeros(1,agents);              
% Put "1" where agents are "alive", then randomize the array
for c=1:Sa(1),                      tempA(c)=1;             end
tempA = RandArray(tempA);           % Randomize array
Av = [tempA ; tempA];               % Initialize vector storing state of A agents     
Bv = ~Av;                           % Initialize vector storing state of B agents 
% Notes on Av, Bv:
% - Markov process, so only previous and current time steps needed --> 2 rows:          
% - Av and Bv are complementary (A+B=agents)
clear c tempA tempB;

%% ABK Simulation

while t*dt <= maxTime
%     dt = 1 / (ko*Sa(t-1));            % Variable time step increment ; OLD
%     dt = exprnd(1/(k*S(t-1)));   % Exponentially-distributed time step increment ; OLD
    	
    for i = 1:agents
        if Av(1,i) == 1                     % if A agent is alive
%             P(t) = 1 - exp(-ko_f * dt);      % P from integrated rate law
            P(t) = ko_f * dt;                % P from differential rate law
            if rand < P(t)                  % **check probability condition**
                Av(2,i) = 0;                 % A agent dies
                Bv(2,i) = 1;                 % A is converted to B
            end
        elseif Bv(1,i) == 1					% if B agent is alive
%        	  P(t) = 1 - exp(-ko_r * dt);      % P from integrated rate law
            P(t) = ko_r * dt;                % P from differential rate law
            if rand < P(t)                  % **check probability condition**
                Bv(2,i) = 0;                 % B agent dies
                Av(2,i) = 1;                 % B is converted to A
            end
        end
    end
    Sa(t+1) = sum(Av(2,:));                  Sb(t+1) = sum(Bv(2,:));
    Av(1,:) = Av(2,:);                       Bv(1,:) = Bv(2,:); 
    t = t + 1;   
end

% Remove unnecessary terminal 0's from arrays
% Sa = Sa(Sa~=0);                             w = size(Sa,2);                 
% Sb = Sb(1:w);                               P = P(1:w); 

disp(['ABK sim terminal A value   = ' num2str(Sa(end))]);      
clear Av Bv;
finaltime = (t-1) * dt;

%% Solve differential equation (was waiting for variable 'finaltime')
[t_sol, y_sol] = ode45(@o1_rev_dif,0:finaltime/100:finaltime+(finaltime/20),[Sa(1) ; Sb(1)]);
disp(['Diff eq terminal A value   = ' num2str(y_sol(end,1))]);
% Note: extended time domain of solution by finaltime/20 to get closer to steady-state values

%% Graph ABM (stochastic) and deterministic results
figure('Name','1st Order Reversible Rx','NumberTitle','off');
time = 0:dt:finaltime;
scatter(time,Sa,3,'.c');                                        hold on;
scatter(time,Sb,3,'.m');
scatter(t_sol,y_sol(:,1),3,'.b');                  
scatter(t_sol,y_sol(:,2),3,'.r');                               hold off;
axis([0 t_sol(end) 0 agents]);          xlabel('time');         ylabel('# Agents');         
legend('A stoc','B stoc','A deter','B deter');

%% Finish
clear r w i j temp;
toc
