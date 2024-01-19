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

% **************************************************************
% Considering heterogeneity in the regulation constant Km_r.
% **************************************************************
% Comparing average trajectories of heterogeneous populations
% (wrt Km_r) using the Common Random Numbers (CRN) procedure.
% Effectively performing sensitivity analysis on the index of
% heterogeneity, psi.
% **************************************************************

% This script is adapted from '../mAct_hetKm_reps.m' 

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear;     close all;     clc;     tic;     

reps = 1000;                       % Number of times to repeat experiment

totalTime = 5000;                  % Simulation time (sec)
dt = 1/50;                         % Constant (fixed) time step increment (sec)
t_steps = totalTime / dt;

agents = 100;                      % Max number of agents for any species (choose even number)
k_b = 0.02;                        % 0th order R synthesis rate, UPregulated by Ep
k_d = 0.075;                       % 1st order degradation rate (units: 1/sec)
k_s = 0.05;                        % 0th order R synthesis rate UPregulated by S
k_f = 1;                           % basal forward rate (1st order)
k_r = 0.05;                        % basal reverse rate (1st order)
Km_f = 5;                          % MICROSCOPIC Michaelis-Menten constant for forward rx
Km_r_mean = 10;                    % MICROSCOPIC Michaelis-Menten constant for reverse rx

S = 15;                            % Assume number of S molecules/agents is NOT changing

% Initial population sizes
Ri = 0;             Epi = 0;        Ei = agents - Epi;					

% Preallocate memory for simulation variables
ABK_R_ss = zeros(reps,2);
R_all = zeros(reps, t_steps+1,2);
Ep_all = zeros(reps, t_steps+1,2);


for n = 1:reps
  fprintf('.');
  
  for m=1:2     % Run simulation for homogeneous and heterogeneous model variants
    rng(n);     % Ensure each run uses the same sequence of random numbers

    % Agent-specific Km_r values for species E. 
    if m==1                            
        % Two subpopulations of equal size with Km_r=1,19 (mean=10, sdev=9.0453)
        Km_r = [ones(1,agents/2) , 19*ones(1,agents/2)];
        het1 = het_Stats(Km_r);     % structure with calculated measures of heterogeneity
    else                              
        % Two subpopulations of (marginally) unequal size with Km_r=1,19 (mean=9.82, sdev=9.0435)
        Km_r = [ones(1,agents/2 +1) , 19*ones(1,agents/2 -1)]; 
        het2 = het_Stats(Km_r);     % structure with calculated measures of heterogeneity
    end
     
    Tr = zeros(1,t_steps);
    Tep = zeros(1,t_steps);             Te = zeros(1,t_steps);
    P_f = zeros(1,t_steps);             P_r = zeros(agents,t_steps);
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
            P_r(tempE(j),t) = k_r * Tr(t) / ( Km_r(tempE(j)) + Te(t) ) * dt;   % wrt each E molecule
            if rand < P_r(tempE(j),t)
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

    R_all(n,:,m) = Tr;                         Ep_all(n,:,m) = Tep;

    ABK_R_ss(n,m) = mean(Tr(750/dt:end));        % Avg R pop size from t=750-->totalTime
  end               % end 'for m'  
end                 % end 'for n'


finaltime = (t-1) * dt;
time = 0:dt:finaltime;

%% Calculate averages
avgR_het1 = mean(R_all(:,:,1));              avgEp_het1 = mean(Ep_all(:,:,1));
avgR_het2 = mean(R_all(:,:,2));              avgEp_het2 = mean(Ep_all(:,:,2));

sdevR_het1 = std(R_all(:,:,1));              sdevEp_het1 = std(Ep_all(:,:,1));
sdevR_het2 = std(R_all(:,:,2));              sdevEp_het2 = std(Ep_all(:,:,2));

%% Start here when loading saved data file ('matlab.mat')
%% Plot selected trajectories: HET1 vs HET2
fig0 = figure('Name','Mutual Activation Time Course','NumberTitle','off');
set(fig0,'Position',[501 1 500 450]);                                   hold on;
set(fig0,'PaperPosition',[0 0 6 5],'PaperUnits','inches');

trial = 1;                      % specify ABK simulation trial #   
% (Interesting switch behavior [Run1]: # 3,4,6,7,8,9,10,12,14,15,16,18,19,27,34,44,47,49,50
%           for Run1, totalTime = 10000, dt = 0.01 sec;
%           for Run2, totalTime = 5000 , dt = 0.02 sec. )

gap = 20;

p01(1) = plot(time(1:gap:end),Ep_all(trial,1:gap:end,1),':','MarkerSize',2,...
    'DisplayName','$\textrm{HET1 } N_{Ep}(t)$');                                          
p01(2) = plot(time(1:gap:end),R_all(trial,1:gap:end,1),...
    'DisplayName','$\textrm{HET1 } N_{R}(t)$');

p02(1) = plot(time(1:gap:end),Ep_all(trial,1:gap:end,2),':','MarkerSize',2,...
    'DisplayName','$\textrm{HET2 } N_{Ep}(t)$');                                         
p02(2) = plot(time(1:gap:end),R_all(trial,1:gap:end,2),...
    'DisplayName','$\textrm{HET2 } N_{R}(t)$');

axis([0 time(end) 0 110]);                              
% axis tight;
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('$t \; \textrm{(sec)}$','Interpreter','LaTeX','FontSize',11);                 
ylabel('$N(t)$','Interpreter','LaTeX','FontSize',11);                             hold off;

leg0 = legend([p01 , p02]);
set(leg0,'Interpreter','LateX','FontSize',10,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');

%% Plot SENSITIVITY of selected (trial) trajectory wrt index of heterogeneity, psi
fig1 = figure('Name','Sensitivity of selected trajectory','NumberTitle','off');
set(fig1,'Position',[501 1 500 450]);                                   % hold on;
set(fig1,'PaperPosition',[0 0 6 5],'PaperUnits','inches');

gap = 20;
sens_R = ( R_all(trial,1:gap:end,1) - R_all(trial,1:gap:end,2) ) / ( het1.psi - het2.psi );
sens_Ep = ( Ep_all(trial,1:gap:end,1) - Ep_all(trial,1:gap:end,2) ) / ( het1.psi - het2.psi );

p1(1) = plot(time(1:gap:end),sens_R * 10^-5,'DisplayName','$\textrm{Species } R$');
p1(2) = plot(time(1:gap:end),sens_Ep * 10^-5,'DisplayName','$\textrm{Species } Ep$');

xlim([0 time(end)]);                              
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('$t \; \textrm{(sec)}$','Interpreter','LaTeX','FontSize',11);                 
ylabel('$sens \{ N_R(t;\,\psi_i),\epsilon \} \, \times 10^{-5}$',...
    'Interpreter','LaTeX','FontSize',11);                             

leg1 = legend(p1);          % leads to error. why?
set(leg1,'Interpreter','LateX','FontSize',10,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');

%% Plot average trajectories: HET1 vs HET2
fig2 = figure('Name','Mutual Activation Average Time Course','NumberTitle','off');
set(fig2,'Position',[501 1 500 450]);                                   hold on;
set(fig2,'PaperPosition',[0 0 6 5],'PaperUnits','inches');

gap = 20;

p1(1) = plot(time(1:gap:end),avgEp_het1(1:gap:end),'LineWidth',1,...
    'DisplayName','$\textrm{HET1 } < \! N_{Ep}(t) \! >$'); 
p1(2) = plot(time(1:gap:end),avgR_het1(1:gap:end),'LineWidth',1,...
    'DisplayName','$\textrm{HET1 } < \! N_{R}(t) \! >$'); 

p2(1) = plot(time(1:gap:end),avgEp_het2(1:gap:end),'LineWidth',1,...
    'DisplayName','$\textrm{HET2 } < \! N_{Ep}(t) \! >$');
p2(2) = plot(time(1:gap:end),avgR_het2(1:gap:end),'LineWidth',1,...
    'DisplayName','$\textrm{HET2 } < \! N_{R}(t) \! >$');

xlim([0 time(end)]);                              
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('$t \; \textrm{(sec)}$','Interpreter','LaTeX','FontSize',11);                 
ylabel('$N(t)$','Interpreter','LaTeX','FontSize',11);                             

leg2 = legend([p1(2) , p2(2) , p1(1) , p2(1)]);
set(leg2,'Interpreter','LateX','FontSize',10,...
    'EdgeColor',[0.95 0.95 0.95],'Location','SouthEast');

%% Plot SDEV of average trajectories: HET1 vs HET2
fig2a = figure('Name','Mutual Activation SDev Time Course','NumberTitle','off');
set(fig2a,'Position',[501 1 500 450]);                                   hold on;
set(fig2a,'PaperPosition',[0 0 6 5],'PaperUnits','inches');

gap = 20;

p1a(1) = plot(time(1:gap:end),sdevEp_het1(1:gap:end),'LineWidth',1,...
    'DisplayName','$\textrm{HET1 } SDev\left\{ < \! N_{Ep}(t) \! > \right\}$',...
    'LineStyle',':'); 
p1a(2) = plot(time(1:gap:end),sdevR_het1(1:gap:end),'LineWidth',1,...
    'DisplayName','$\textrm{HET1 } SDev\left\{ < \! N_{R}(t) \! > \right\}$',...
    'LineStyle',':'); 

p2a(1) = plot(time(1:gap:end),sdevEp_het2(1:gap:end),'LineWidth',1,...
    'DisplayName','$\textrm{HET2 } SDev\left\{ < \! N_{Ep}(t) \! > \right\}$',...
    'LineStyle',':');
p2a(2) = plot(time(1:gap:end),sdevR_het2(1:gap:end),'LineWidth',1,...
    'DisplayName','$\textrm{HET2 } SDev\left\{ < \! N_{R}(t) \! > \right\}$',...
    'LineStyle',':');

xlim([0 time(end)]);                              
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('$t \; \textrm{(sec)}$','Interpreter','LaTeX','FontSize',11);                 
ylabel('$SDev\left\{ < \! N(t) \! > \right\}$','Interpreter','LaTeX','FontSize',11);                             

leg2a = legend([p1a(2) , p2a(2) , p1a(1) , p2a(1)]);
set(leg2a,'Interpreter','LateX','FontSize',10,...
    'EdgeColor',[0.95 0.95 0.95],'Location','SouthEast');


%% Plot SENSITIVITY of average trajectories wrt index of heterogeneity, psi
fig3 = figure('Name','Sensitivity of average trajectories','NumberTitle','off');
set(fig3,'Position',[501 1 500 450]);                                           hold on;
set(fig3,'PaperPosition',[0 0 6 5],'PaperUnits','inches');

gap = 20;

sens_avgR = ( avgR_het1(1:gap:end) - avgR_het2(1:gap:end) ) / ( het1.psi - het2.psi );
sens_avgEp = ( avgEp_het1(1:gap:end) - avgEp_het2(1:gap:end) ) / ( het1.psi - het2.psi );

p3(1) = plot(time(1:gap:end),sens_avgR * 10^-4,...
    'DisplayName','$\textrm{Species } R$',...
    'Color',[0.467, 0.675, 0.188]);
p3(2) = plot(time(1:gap:end),sens_avgEp * 10^-4,...
    'DisplayName','$\textrm{Species } Ep$',...
    'Color',[0.494, 0.184, 0.557]);

plot([0 time(end)],[0 0],'--','LineWidth',2,'Color','k');

xlim([0 time(end)]);  
ylim([-1.15 1.15]);

set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('$t \; \textrm{(sec)}$','Interpreter','LaTeX','FontSize',11);                 
ylabel('$sens \{ < \! N(t;\,\psi_i) \! > ,\, \epsilon \} \, \times 10^{-4}$',...
    'Interpreter','LaTeX','FontSize',11);           

leg3 = legend(p3);
set(leg3,'Interpreter','LateX','FontSize',10,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');

%% Plot SENSITIVITY of SDEV of average trajectories wrt index of heterogeneity, psi
fig3a = figure('Name','Sensitivity of average trajectories','NumberTitle','off');
set(fig3a,'Position',[501 1 500 450]);                                           hold on;
set(fig3a,'PaperPosition',[0 0 6 5],'PaperUnits','inches');

gap = 20;

sens_sdevR = ( sdevR_het1(1:gap:end) - sdevR_het2(1:gap:end) ) / ( het1.psi - het2.psi );
sens_sdevEp = ( sdevEp_het1(1:gap:end) - sdevEp_het2(1:gap:end) ) / ( het1.psi - het2.psi );

p3a(1) = plot(time(1:gap:end),sens_sdevR * 10^-4,...
    'DisplayName','$\textrm{Species } R$',...
    'LineStyle',':',...
    'Color',[0.467, 0.675, 0.188]);
p3a(2) = plot(time(1:gap:end),sens_sdevEp * 10^-4,...
    'DisplayName','$\textrm{Species } Ep$',...
    'LineStyle',':',...
    'Color',[0.494, 0.184, 0.557]);

plot([0 time(end)],[0 0],'--','LineWidth',2,'Color','k');

xlim([0 time(end)]);  
ylim([-1.3 1.3]);

set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
yticks([-1.2:0.2:1.2]);

xlabel('$t \; \textrm{(sec)}$','Interpreter','LaTeX','FontSize',11);                 
ylabel('$sens \{ SDev \{ < \! N(t;\,\psi_i) \! > \} ,\, \epsilon \} \, \times 10^{-4}$',...
    'Interpreter','LaTeX','FontSize',11);           

leg3a = legend(p3a);
set(leg3a,'Interpreter','LateX','FontSize',10,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');

%% Detect and plot the system's switch state for all trajectories.

% Fixed points calculated for a system with homogeneous populations.
% The following values apply to the parameters used in the report.
Ep_ss_ON = 86.69;       % Calculated elsewhere. Change accordingly.
Ep_ss_OFF = 5.35;       % Calculated elsewhere. Change accordingly.

ON = detect_SwitchState(Ep_all,time, 50, [Ep_ss_ON, Ep_ss_OFF], ["HET1","HET2"]);

%% % Saving all trajectories takes up too much space. Keep just the ones 
% that have been manually screened and selected for interesting behavior.
% Use this as necessary.

% keep_Trajs = [3,4,6,7,8,9,10,12,14,15,16,18,19,27,34,44,47,49,50];
% 
% for i=size(R_all,1):-1:1
%     if find(keep_Trajs==i) == 1
%         continue;
%     else
%         R_all(i,:,:) = [];
%         Ep_all(i,:,:) = [];
%     end
% end

%% Finish - Notes
clear h i j k m n t Tr Tep Te R Ep E temp tempR tempEp tempE fig* leg* p1 p2 p3;
toc
