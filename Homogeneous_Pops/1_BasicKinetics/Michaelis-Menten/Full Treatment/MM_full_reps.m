% E + A <--> EA --> E + B
% Simulating Michaelis-Menten kinetics 
% (Full treatment: all three rxs are explicitly modeled)
% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear;            tic;                  clc; 
rng(0);
pb = waitbar(0,'0');
%% Declare variables and functions
global k_f k_r k_cat;

maxTime = 300;                          % Maximum Simulation time (sec)
dt = 1/100;                             % Constant (fixed) time step increment (sec)
steps = maxTime / dt;

reps = 500;
% Preallocate matrices storing results of all simulation repetitions
Ta = zeros(reps,steps);                 Tb = zeros(reps,steps);
Te = zeros(reps,steps);                 Tea = zeros(reps,steps);

agents = 20;                disp(['Agents = ' num2str(agents)]);
Ao = agents;                Bo = 0;             EAo = 0;
Et = 5;                                 % Total number of enzyme molecules

k_f = 0.01;                             % MICROscopic forward rate (2nd order)
k_r = 0.05;                              % reverse rx rate (1st order)
k_cat = 0.05;                            % catalytic rx rate EA --> E + A (1st order)
Km = (k_r + k_cat) / k_f;               % Michaelis-Menten constant (Microscopic);
%% Start reps
for n=1:reps
%     fprintf(1,'.');
    progress = n/reps;
    waitbar(progress,pb,sprintf('%.0f%% %',progress*100));

    % Initialize - Preallocate memory for variables; ** Initial conditions **
    Pa = zeros(1,steps);              % For storing probability value for rx of A at each time step
    Pea = zeros(1,steps);             % For storing probability value for rx of EA at each time step
    
    Sa = zeros(1,steps);              % Sum of A agents in each time step	
    Sb = zeros(1,steps);              % Sum of B agents in each time step
    Se = zeros(1,steps);              % Sum of E (enzyme) agents in each time step
    Sea = zeros(1,steps);             % Sum of EA (EA complex) agents in each time step

    % ******** Initial conditions - Number of A, B Agents ********
    Sa(1) = Ao;                         Sb(1) = Bo;
    Se(1) = Et;                         Sea(1) = EAo;
    % ************************************************************

    Av = ones(2,agents);                % Initialize vector storing state of A agents     
    Bv = zeros(1,agents);               % Initialize vector storing state of B agents 
    Ev = ones(2,Et);                    % Initialize vector storing state of E (enzyme) agents
    EAv = ~Ev;                          % Initialize vector storing state of EA (EA complex) agents
    % Notes on Av, Bv:
    % - Markov process, so only previous and current time steps needed --> 2 rows
    % - B molecules are not reactants in a rx so only one row is need to keep track of them.
    % - Ev and EAv are complementary (E + EA = Et)

    % ABM Simulation
    time = zeros(1,steps);
    
    for t=2:steps                % while Sa(t) > 1 && t*dt <= maxTime

        for i = 1:agents
            if Av(1,i) == 1                      % if A agent is alive
                Pa(t-1) = k_f * Se(t-1) * dt;        % P from differential rate law
                if rand < Pa(t-1)                  % **check probability condition**
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
        Pea(t-1) = h * dt;                         % P from differential rate law
        nk_r = k_r / h;                          % Normalized k_r value
        nk_cat = k_cat / h;                      % Normalized k_cat value

        for j = 1:Et
            if EAv(1,j) == 1                     % if EA agent is alive
                r = rand;
                if r < nk_r * Pea(t-1)             % **Reverse Rx probability condition**
                    EAv(2,j) = 0;                % EA agent "dies" (E + A <-- EA)
                    Ev(2,j) = 1;                 % EA is converted to E
                    temp1 = find(Av(1,:)==0);    % Randomly choose A agent which is "born"
                    y = ceil(rand * size(temp1,2));
                    Av(2,temp1(y)) = 1;         
                elseif r >= nk_r * Pea(t-1) && r < Pea(t-1)  % **Forward (cat) Rx probability condition**
                    EAv(2,j) = 0;                % EA agent "dies" (EA --> E + B)
                    Ev(2,j) = 1;                 % EA is converted to E
                    z = find(Bv==0,1);           % Find first occurence of B agent who is currently "dead"
                    Bv(z)= 1;                    % B agent is "born"
                end
            end
        end

        Sa(t) = sum(Av(2,:));                  Sb(t) = sum(Bv);
        Se(t) = sum(Ev(2,:));                  Sea(t) = sum(EAv(2,:));

        Av(1,:) = Av(2,:);                       
        Ev(1,:) = Ev(2,:);                       EAv(1,:) = EAv(2,:);
        
        time(t) = time(t-1) + dt;
    end             % end "for t" loop
    Ta(n,:) = Sa;                       Tb(n,:) = Sb;
    Te(n,:) = Se;                       Tea(n,:) = Sea;    
end                 % end "n=1:reps" for loop

avgA = mean(Ta);                     sdevA = std(Ta);
avgB = mean(Tb);                     sdevB = std(Tb);
avgE = mean(Te);                     sdevE = std(Te);
avgEA = mean(Tea);                   sdevEA = std(Tea);

clear Av Bv Ev EAv nk_r nk_cat;

%% Solve differential equation
[t_sol, y_sol] = ode45(@mm_full_dif,time,[Ao ; Bo; Et; EAo]);  % Same time sampling as in simulation
% Or, uncomment the following to try a different time sampling
% [t_sol, y_sol] = ode45(@mm_full_dif,0:maxTime/100:maxTime,[Ao ; Bo; Et; EAo]);

%% Find Coefficient of Determination, R^2
% This works for any curve fit, linear or nonlinear.
SST(1) = sum((avgA' - mean(avgA)).^2);     % Total sum of squares for simulation data (for A)
SSR(1) = sum((avgA' - y_sol(:,1)).^2);     % sum of square residuals (sim data vs DE predictions)

SST(2) = sum((avgB' - mean(avgB)).^2);     % Total sum of squares for simulation data (for B)
SSR(2) = sum((avgB' - y_sol(:,2)).^2);     % sum of square residuals (sim data vs DE predictions)

SST(3) = sum((avgE' - mean(avgE)).^2);     % Total sum of squares for simulation data (for E)
SSR(3) = sum((avgE' - y_sol(:,3)).^2);     % sum of square residuals (sim data vs DE predictions)

SST(4) = sum((avgEA' - mean(avgEA)).^2);   % Total sum of squares for simulation data (for EA)
SSR(4) = sum((avgEA' - y_sol(:,4)).^2);    % sum of square residuals (sim data vs DE predictions)

Rsq = 1 - SSR./SST                         % Definition of R^2
R = sqrt(Rsq);                             % Correlation Coefficient R

% Reminder: % E + A <--> EA --> E + B

%% Plot Time Courses
figure1 = figure('Name','Michaelis-Menten (full) TC','NumberTitle','off');
set(figure1,'Position',[1 1 500 450]);
p1 = plot(time,avgA,'b','MarkerSize',3,'DisplayName','<N_A(t)>_{sim}');                 hold on;
% p1_dev1 = plot(time,avgA+sdevA,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p1_dev0 = plot(time,avgA-sdevA,'LineStyle','--','Color',[0.8 0.8 0.8]);
p1s = plot(time,Ta(end-5,:),'m','MarkerSize',3,'DisplayName','Sample Trajectory');

p2 = plot(time,avgB,'r','MarkerSize',3,'DisplayName','<N_X(t)>_{sim}');                 
p2_dev1 = plot(time,avgB+sdevB,'LineStyle','--','Color',[0.8 0.8 0.8]);
p2_dev0 = plot(time,avgB-sdevB,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p2s = plot(time,Tb(end,:),'m','MarkerSize',3,'DisplayName','Sample Trajectory');

p3 = plot(time,avgE,'g','MarkerSize',3,'DisplayName','<N_E(t)>_{sim}');                 
% p3_dev1 = plot(time,avgE+sdevE,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p3_dev0 = plot(time,avgE-sdevE,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p3s = plot(time,Te(end-1,:),'m','MarkerSize',3,'DisplayName','Sample Trajectory');

p4 = plot(time,avgEA,'k','MarkerSize',3,'DisplayName','<N_{EA}(t)>_{sim}');                 
p4_dev1 = plot(time,avgEA+sdevEA,'LineStyle','--','Color',[0.8 0.8 0.8]);
p4_dev0 = plot(time,avgEA-sdevEA,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p4s = plot(time,Tea(end,:),'m','MarkerSize',3,'DisplayName','Sample Trajectory');

p1d = plot(t_sol,y_sol(:,1),':b','DisplayName','DE');        
p2d = plot(t_sol,y_sol(:,2),':r');
p3d = plot(t_sol,y_sol(:,3),':g');
p4d = plot(t_sol,y_sol(:,4),':k');

xlabel('t (sec)');                 ylabel('N(t)');                  hold off;   
axis([0 maxTime 0 agents]);          % axis tight;                                 
% title(,'FontName','Times New Roman','FontSize',11)
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
leg = legend([p1 p2 p3 p4]);
% set(leg,'OuterPosition',[0.673 0.374 0.182 0.123]);
set(leg,'Location','East');
set(leg,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95]);
% set(p??,'Visible','off');
%% Finish
clear h i j n r t w x y z temp temp1 p* figure1 leg;
toc
