% E + A <--> EA --> E + B
% Simulating Michaelis-Menten kinetics 
% (Full treatment: all three rxs are explicitly modeled)
% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;            tic;
clc; 
rng(0);
%% Declare variables and functions
global k_f k_r k_cat;

maxTime = 10000;                        % Maximum Simulation time (sec)
dt = 1/100;                             % Constant (fixed) time step increment (sec)
t_steps = maxTime / dt;

agents = 50;                disp(['Agents = ' num2str(agents)]);
Et = 5;                                 % Total number of enzyme molecules

k_f = 0.005;                             % MICROscopic forward rate (2nd order)
k_r = 0.05;                              % reverse rx rate (1st order)
k_cat = 0.05;                            % catalytic rx rate EA --> E + A (1st order)
Km = (k_r + k_cat) / k_f;               % Michaelis-Menten constant (Microscopic);
%% Initialize - Preallocate memory for variables; ** Initial conditions **
t = 1;                              % Time counter variable

Pa = zeros(1,t_steps);              % For storing probability value for rx of A at each time step
Pea = zeros(1,t_steps);             % For storing probability value for rx of EA at each time step
Sa = zeros(1,t_steps);              % Sum of A agents in each time step	

Sb = zeros(1,t_steps);              % Sum of B agents in each time step
Se = zeros(1,t_steps);              % Sum of E (enzyme) agents in each time step
Sea = zeros(1,t_steps);             % Sum of EA (EA complex) agents in each time step

% ******** Initial conditions - Number of A, B Agents ********
Sa(1) = agents;                		Sb(1) = 0;
Se(1) = Et;                         Sea(1) = 0;
% ************************************************************

Av = ones(2,agents);                % Initialize vector storing state of A agents     
Bv = zeros(1,agents);               % Initialize vector storing state of B agents 
Ev = ones(2,Et);                    % Initialize vector storing state of E (enzyme) agents
EAv = ~Ev;                          % Initialize vector storing state of EA (EA complex) agents
% Notes on Av, Bv:
% - Markov process, so only previous and current time steps needed --> 2 rows
% - B molecules are not reactants in a rx so only one row is need to keep track of them.
% - Ev and EAv are complementary (E + EA = Et)

%% ABM Simulation

while Sa(t) > 1 && t*dt <= maxTime
    	
    for i = 1:agents
        if Av(1,i) == 1                      % if A agent is alive
            Pa(t) = k_f * Se(t) * dt;        % P from differential rate law
            if rand < Pa(t)                  % **check probability condition**
                temp = find(Ev(1,:)==1);     % Randomly choose E agent which forms EA complex
                if isempty(temp) == 1
                    break;
                else
                Av(2,i) = 0;                 % A is converted to EA
                x = ceil(rand * size(temp,2));
                Ev(2,temp(x)) = 0;           % E (enzyme) agent becomes occupied
                EAv(2,temp(x)) = 1;          % EA complex forms              
                end
            end
        end
    end
    
    h = k_r + k_cat;                         % total rate of conversion of EA
    Pea(t) = h * dt;                         % P from differential rate law
    nk_r = k_r / h;                          % Normalized k_r value
    nk_cat = k_cat / h;                      % Normalized k_cat value
       
    for j = 1:Et
        if EAv(1,j) == 1                     % if EA agent is alive
            r = rand;
            if r < nk_r * Pea(t)             % **Reverse Rx probability condition**
                EAv(2,j) = 0;                % EA agent "dies" (E + A <-- EA)
                Ev(2,j) = 1;                 % EA is converted to E
                temp1 = find(Av(1,:)==0);    % Randomly choose A agent which is "born"
                y = ceil(rand * size(temp1,2));
                Av(2,temp1(y)) = 1;         
            elseif r >= nk_r * Pea(t) && r < Pea(t)  % **Forward (cat) Rx probability condition**
                EAv(2,j) = 0;                % EA agent "dies" (EA --> E + B)
                Ev(2,j) = 1;                 % EA is converted to E
                z = find(Bv==0,1);           % Find first occurence of B agent who is currently "dead"
                Bv(z)= 1;                    % B agent is "born"
            end
        end
    end
            
    Sa(t+1) = sum(Av(2,:));                  Sb(t+1) = sum(Bv);
    Se(t+1) = sum(Ev(2,:));                  Sea(t+1) = sum(EAv(2,:));
    
    Av(1,:) = Av(2,:);                       
    Ev(1,:) = Ev(2,:);                       EAv(1,:) = EAv(2,:);
    t = t + 1;   
end

% Remove unnecessary terminal 0's from arrays
Sa = Sa(Sa~=0);                             w = size(Sa,2);                 
Sb = Sb(1:w);                               Se = Se(1:w);               Sea = Sea(1:w);                             
Pa = Pa(1:w);                               Pea = Pea(1:w);          
clear Av Bv Ev EAv;
finaltime = (t-1) * dt;

%% Solve differential equation
[t_sol, y_sol] = ode45(@mm_full_dif,0:finaltime/100:finaltime,[Sa(1) ; Sb(1); Se(1); Sea(1)]);

%% Graph ABK (stochastic) time course
figure('Name','Michaelis-Menten Rx Time course','NumberTitle','off');               hold on;
time = 0:dt:finaltime;
plot(time,Sa,'b');                           plot(time,Sb,'r');
plot(time,Se,'g');                           plot(time,Sea,'m');
plot(t_sol,y_sol(:,1),':b');               plot(t_sol,y_sol(:,2),':r');
plot(t_sol,y_sol(:,3),':g');               plot(t_sol,y_sol(:,4),':m'); 
axis([0 finaltime 0 agents]);          
xlabel('time');          ylabel('# Agents');         
legend('A','B','E','EA','Location','Best');                                                           hold off;

%% Finish
clear h i j r t w x y z temp temp1;
toc

%% Notes
% - Note that agreement with differential equation is best during the initial reaction time period
% because that's when the EA steady state approximation is valid, as predicted.
