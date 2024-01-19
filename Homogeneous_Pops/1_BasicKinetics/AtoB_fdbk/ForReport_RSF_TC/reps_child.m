%  A --> B  ,     with feedback from B
% Simulating 1st order kinetics. Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

% clear;  clc;     tic;

global ko n K;
ko = 0.1;            % basal rate
% Feedback parameters
% alpha = 0;           % degree of activation: >1 activator, <1 repressor, =1 no regulation
n = 2;               % Hill Coefficient: measure of cooperativity
K = 10;              % Agents needed for half-maximal Actvation/Repression

agents = 20;

t_max = 50;             % Simulation time (sec)
dt = 1/100;                % Constant (fixed) time step increment (sec)
steps = t_max / dt;

reps = 500;

P = zeros(reps,steps);
Sa = zeros(reps,steps);           Sb = zeros(reps,steps);

for z=1:reps

    % Initial number of reactant molecules
    Ao = agents;                Bo = 0;        
    % Arrays for tracking reactant agent status
    A = ones(1,Ao);             B = zeros(1,Ao);               
    % Initialize time-dependent sum of molecule numbers
    At = zeros(1,steps);        Bt = zeros(1,steps);   
    At(1) = sum(A);             Bt(1) = sum(B);               

    time = zeros(1,steps);
    
    for t = 2:steps
    %     dt = 1 / (ko*Sa(t-1));            % Variable time step increment
    %     dt = exprnd(1/(k*S(t-1)));   % Exponentially-distributed time step increment

        F =  (K^n + alpha * Bt(t-1)^n) / (K^n + Bt(t-1)^n);   % RSF
        k = ko * F;

%         P(z,t-1) = 1 - exp(-k*dt);      % P from integrated rate law
        P(z,t-1) = k * dt;                  % P from differential rate law

        tempA = find(A==1);
         
        for i = 1:size(tempA,2)
            if rand < P(z,t-1)              % **check probability condition**
                A(tempA(i)) = 0;            % A agent dies
                B(tempA(i)) = 1;            % B is born
            end
        end   
        At(t) = sum(A);                 Bt(t) = sum(B);
        time(t) = time(t-1) + dt;
    end         % "for t"
    Sa(z,:) = At;                   Sb(z,:) = Bt;
end             % "for z" (reps)

avgA = mean(Sa);                     sdevA = std(Sa);
avgB = mean(Sb);                     sdevB = std(Sb);
%% Solve ODE for mixed order kinetics
[t_sol,y_sol] = ode45(@o1f_simple_dif,0:t_max/500:t_max,[Ao ; Bo]);
%% Plot Time Courses
% figure1 = figure('Name','o1 rx w/ Feedback ','NumberTitle','off');
% set(figure1,'Position',[1 1 500 450]);
% p1 = plot(time,avgA,'b','MarkerSize',3,'DisplayName','<N_A(t)>_{sim}');                 hold on;
% p1_dev1 = plot(time,avgA+sdevA,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p1_dev0 = plot(time,avgA-sdevA,'LineStyle','--','Color',[0.8 0.8 0.8]);
% % p1s = plot(time,Sa(end,:),'m','MarkerSize',3,'DisplayName','Sample Trajectory');
% 
% p2 = plot(time,avgB,'r','MarkerSize',3,'DisplayName','<N_B(t)>_{sim}');                 
% p2_dev1 = plot(time,avgB+sdevB,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p2_dev0 = plot(time,avgB-sdevB,'LineStyle','--','Color',[0.8 0.8 0.8]);
% % p2s = plot(time,Sb(end,:),'m','MarkerSize',3,'DisplayName','Sample Trajectory');
% 
% p1d = plot(t_sol,y_sol(:,1),':b','DisplayName','DE');        
% p2d = plot(t_sol,y_sol(:,2),':r');
% 
% xlabel('t (sec)');                 ylabel('N(t)');                  hold off;   
% axis tight;                         % axis([0 t_max 0 agents]); 
% set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
% leg = legend([p1 p2]);
% set(leg,'Location','NorthEast');
% set(leg,'FontName','Times New Roman','FontSize',9,...
%     'EdgeColor',[0.95 0.95 0.95]);

%% Finish
clear tempA t w i z;
% toc
