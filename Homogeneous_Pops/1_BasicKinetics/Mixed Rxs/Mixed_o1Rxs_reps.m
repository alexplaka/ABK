%  A     --> C (or X)             Rate constant k1  (1st order)
%  A     --> D (or Y)             Rate constant k2  (1st order)
% Simulating concurrent 1st order reactions using my agent-based
% algorithm. 

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear; clc;  tic
rng(0);

global k1 k2;
k1 = 0.5;                % k1: 1st order rate constant [units: 1/sec]
k2 = 0.2;                % k2: 1st order rate constant [units: 1/sec]

t_max = 5;                             % in seconds
dt = 1/100;                             % Fixed time step increment
steps = t_max / dt;

reps = 500;
Sa = zeros(reps,steps);                 
Sc = zeros(reps,steps);                 Sd = zeros(reps,steps);

for n=1:reps

    % Initial number of reactant molecules; Assume B is limiting for 2nd order rx.
    Ao = 100;            
    % Arrays for tracking reactant agent status
    A = ones(1,Ao);                           
    % Initialize time-dependent sum of molecule numbers
    At = zeros(1,steps);
    At(1) = sum(A); 
    Ct = zeros(1,steps);            Dt = zeros(1,steps);    
    Ct(1) = 0;                      Dt(1) = 0;

    time = zeros(1,steps);

    for t=2:steps                   %  while Bt(t-1) > 0 && At(t-1) > 0
        
        deltaC = 0;             deltaD = 0;  

    %     dt = 10 / At(t-1);                        % Variable time step incr. 1
    %     dt = 1/(k1*At(t-1) + k2*At(t-1));         % Variable time step incr. 2
    %     dt = exprnd(1/(k1*At(t-1) + k2*At(t-1))); % Exponentially-distributed time step increment
        
        P1 = k1 * dt;                       % P_dif of 1st order rx A --> C
    %     P1 = 1 - exp(-k1 * dt);           % P_ber of 1st order rx A --> C
        P2 = k2 * dt;                       % P_dif-A of 2nd order rx A + B --> D
%         P2 = 1 - exp(-k2 * dt);           % P_ber-A of 1st order rx A --> D
        
        for i = 1:size(A,2)
            if A(i) == 1      
                r = rand;
                if r < P1                           % 1st order rx A --> C
                    A(i) = 0;           
                    deltaC = deltaC + 1;
                elseif r >= P1 && r < P1+P2         % 1st order rx A --> D
                    A(i) = 0;           
                    deltaD = deltaD + 1;
                end 
            end
        end 
        At(t) = sum(A);                 
        Ct(t) = Ct(t-1) + deltaC;       Dt(t) = Dt(t-1) + deltaD;
        time(t) = time(t-1) + dt;
    end
    Sa(n,:) = At;                       
    Sc(n,:) = Ct;                       Sd(n,:) = Dt;

end         % end 'for n' loop

avgA = mean(Sa);                     sdevA = std(Sa);
avgC = mean(Sc);                     sdevC = std(Sc);
avgD = mean(Sd);                     sdevD = std(Sd);

%% Solve ODE for mixed order kinetics
[t_sol,y_sol] = ode45(@mixed_o1Rxs_dif,time,[At(1); Ct(1); Dt(1)]);  % Same time sampling as in sim
% [t_sol,y_sol] = ode45(@mixed_o1Rxs_dif,0:t_max/100:t_max,[At(1); Ct(1); Dt(1)]); % Coarse time sampling

%% Find Coefficient of Determination, R^2
% This works for any curve fit, linear or nonlinear.
SST(1) = sum((avgA' - mean(avgA)).^2);     % Total sum of squares for simulation data (for A)
SSR(1) = sum((avgA' - y_sol(:,1)).^2);     % sum of square residuals (sim data vs DE predictions)

SST(2) = sum((avgC' - mean(avgC)).^2);     % Total sum of squares for simulation data (for C)
SSR(2) = sum((avgC' - y_sol(:,2)).^2);     % sum of square residuals (sim data vs DE predictions)

SST(3) = sum((avgD' - mean(avgD)).^2);     % Total sum of squares for simulation data (for D)
SSR(3) = sum((avgD' - y_sol(:,3)).^2);     % sum of square residuals (sim data vs DE predictions)

Rsq = 1 - SSR./SST                         % Definition of R^2
R = sqrt(Rsq);                             % Correlation Coefficient R

%% Plot Time Courses
figure1 = figure('Name','Mixed o1 Rx Time course','NumberTitle','off');
set(figure1,'Position',[1 1 500 450]);

p1 = plot(time,avgA,'b','MarkerSize',3,'DisplayName','<N_A(t)>_{sim}');                 hold on;
p1_dev1 = plot(time,avgA+sdevA,'LineStyle','--','Color',[0.8 0.8 0.8]);
p1_dev0 = plot(time,avgA-sdevA,'LineStyle','--','Color',[0.8 0.8 0.8]);
p1s = plot(time,Sa(5,:),'m','MarkerSize',3,'DisplayName','Sample Trajectory');

p3 = plot(time,avgC,'k','MarkerSize',3,'DisplayName','<N_X(t)>_{sim}');                 
p3_dev1 = plot(time,avgC+sdevC,'LineStyle','--','Color',[0.8 0.8 0.8]);
p3_dev0 = plot(time,avgC-sdevC,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p3s = plot(time,Sc(end,:),'Color',[0.5 0.5 0.9],'MarkerSize',3,'DisplayName','Sample Trajectory');

p4 = plot(time,avgD,'g','MarkerSize',3,'DisplayName','<N_Y(t)>_{sim}');                 
p4_dev1 = plot(time,avgD+sdevD,'LineStyle','--','Color',[0.8 0.8 0.8]);
p4_dev0 = plot(time,avgD-sdevD,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p4s = plot(time,Sd(end,:),'Color',[0.1 0 1],'MarkerSize',3,'DisplayName','Sample Trajectory');

p1d = plot(t_sol,y_sol(:,1),':b','DisplayName','DE');        
p3d = plot(t_sol,y_sol(:,2),':k');
p4d = plot(t_sol,y_sol(:,3),':g');

xlabel('t (sec)');                 ylabel('N(t)');                  hold off;                            
axis([0 t_max 0 125]);          % axis tight;                                 
title({['A \rightarrow X  (k_1 = ' num2str(k1) ' sec^{-1})'],...
    ['A \rightarrow Y  (k_2 = ' num2str(k2) ' sec^{-1})']}...
    ,'FontName','Times New Roman','FontSize',11)
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
leg = legend([p1 p3 p4]);
% set(leg,'OuterPosition',[0.673 0.374 0.182 0.123]);
set(leg,'Location','NorthEast');
set(leg,'FontName','Times New Roman','FontSize',9,'EdgeColor',[0.95 0.95 0.95]);
% set(p??,'Visible','off');

%% Finish
clear temp x deltaC deltaD r i p*;
toc

% Result:
% - Works for fixed, variable, and exponentially-distributed time step
% increments. (latter has the smallest computing time).
% - Note that calculated probabilities work best for small dt. For variable
% time step increments some dt's will be large (towards the end of the
% simulation). It may be best to use fixed time step increments while
% calculating such probabilities, despite its computational cost.