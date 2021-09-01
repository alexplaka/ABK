% Simulating Michaelis-Menten kinetics - Simple Implementation
% Overall rx: A --> B
% Actual rxs: E + A <--> EA --> E + B       (Not explicitly modeled)
% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;            tic;                  clc; 
rng(0);
pb = waitbar(0,'0');
%% Declare variables and functions
global k_cat Km E_tot;

maxTime = 300;                          % Maximum Simulation time (sec)
dt = 1/100;                             % Constant (fixed) time step increment (sec)
steps = maxTime / dt;

reps = 500;
% Preallocate matrices storing results of all simulation repetitions
Ta = zeros(reps,steps);                 Tb = zeros(reps,steps);

agents = 20;                disp(['Agents = ' num2str(agents)]);
Ao = agents;                Bo = 0;          
E_tot = 5;                                  % Total number of enzyme molecules
k_cat = 0.05;                               % catalytic rx rate A --> B
Km = 10;                                    % Michaelis-Menten constant (Microscopic);
%% Start reps
for n=1:reps
%     fprintf(1,'.');
    progress = n/reps;
    waitbar(progress,pb,sprintf('%.0f%% %',progress*100));

    % Initialize - Preallocate memory for variables; ** Initial conditions **
    P = zeros(1,steps);               % For storing probability value for rx of A at each time step
    
    Sa = zeros(1,steps);              % Sum of A agents in each time step	

    % ******** Initial conditions - Number of A, B Agents ********
    Sa(1) = Ao;                         Sb(1) = Bo;
    % ************************************************************

    Av = ones(2,agents);                % Initialize vector storing state of A agents     
    % Note on Av: Markov process, so only previous and current time steps needed --> 2 rows

    % ABM Simulation
    time = zeros(1,steps);
    
    for t=2:steps                % while Sa(t) > 1 && t*dt <= maxTime

        for i = 1:agents
            if Av(1,i) == 1                      % if A agent is alive
                
                w = k_cat*E_tot / (Sa(t-1)+Km);
%                 P(t-1) = w * dt;                % P_dif from MM kinetics rate law
                P(t-1) = 1 - exp(-w*dt);        % P_ber
                
                if rand < P(t-1)                  % **check probability condition**
                    
                    Av(2,i) = 0;                 % A is converted to EA
                             
                end
            end
        end

        Sa(t) = sum(Av(2,:));
        Av(1,:) = Av(2,:);                           
        time(t) = time(t-1) + dt;
    end             % end "for t" loop
    Ta(n,:) = Sa;                       Tb(n,:) = agents - Sa;
end                 % end "n=1:reps" for loop

avgA = mean(Ta);                     sdevA = std(Ta);
avgB = mean(Tb);                     sdevB = std(Tb);

%% Solve differential equation
[t_sol, y_sol] = ode45(@mm_sim_dif,time,[Ao ; Bo]);   % Same time sampling as in simulation
% Or, uncomment the following to try a different time sampling
% [t_sol, y_sol] = ode45(@mm_sim_dif,0:maxTime/100:maxTime,[Ao ; Bo]);

%% Find Coefficient of Determination, R^2
% This works for any curve fit, linear or nonlinear.
SST(1) = sum((avgA' - mean(avgA)).^2);     % Total sum of squares for simulation data (for A)
SSR(1) = sum((avgA' - y_sol(:,1)).^2);     % sum of square residuals (sim data vs DE predictions)

SST(2) = sum((avgB' - mean(avgB)).^2);     % Total sum of squares for simulation data (for B)
SSR(2) = sum((avgB' - y_sol(:,2)).^2);     % sum of square residuals (sim data vs DE predictions)

Rsq = 1 - SSR./SST                         % Definition of R^2
R = sqrt(Rsq);                             % Correlation Coefficient R

%% Plot Time Courses
figure1 = figure('Name','Michaelis-Menten (simple) TC','NumberTitle','off');
set(figure1,'Position',[1 1 500 450]);

p1 = plot(time,avgA,'b','MarkerSize',3,'DisplayName','<N_A(t)>_{sim}');                 hold on;
% p1_dev1 = plot(time,avgA+sdevA,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p1_dev0 = plot(time,avgA-sdevA,'LineStyle','--','Color',[0.8 0.8 0.8]);
p1s = plot(time,Ta(end-4,:),'m','MarkerSize',3,'DisplayName','Sample Trajectory');

p2 = plot(time,avgB,'r','MarkerSize',3,'DisplayName','<N_X(t)>_{sim}');                 
p2_dev1 = plot(time,avgB+sdevB,'LineStyle','--','Color',[0.8 0.8 0.8]);
p2_dev0 = plot(time,avgB-sdevB,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p2s = plot(time,Tb(end,:),'m','MarkerSize',3,'DisplayName','Sample Trajectory');

p1d = plot(t_sol,y_sol(:,1),':b','DisplayName','DE');        
p2d = plot(t_sol,y_sol(:,2),':r');

xlabel('t (sec)');                 ylabel('N(t)');                  hold off;   
axis([0 maxTime 0 agents]);          % axis tight;                                 
% title(,'FontName','Times New Roman','FontSize',11)
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
leg = legend([p1 p2]);
% set(leg,'OuterPosition',[0.673 0.374 0.182 0.123]);
set(leg,'Location','East');
set(leg,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95]);
% set(p??,'Visible','off');
%% Finish
close(pb);
clear i t Av pb n Sa Sb p* figure1 leg n;
toc
