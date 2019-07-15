%  A + B     --> X            Rate constant k1  (2nd order)
%  A + C + D --> Y            Rate constant k2  (3rd order)
% Simulating concurrent 2nd and 3rd order kinetics using my agent-based
% algorithm.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear; clc;  tic

pb = waitbar(0,'0');

global k1 k2;
k1 = 0.010;              % k1: 2nd order rate constant [units: 1/sec]
k2 = 0.001;              % k2: 3rd order *microscopic* rate constant [units: 1/sec]

t_max = 1;                              % in seconds
dt = 1/100;                             % Fixed time step increment
steps = t_max / dt;

reps = 500;
Sa = zeros(reps,steps);     Sb = zeros(reps,steps);     
Sc = zeros(reps,steps);     Sd = zeros(reps,steps);
Sx = zeros(reps,steps);     Sy = zeros(reps,steps);

for n=1:reps
%     fprintf(1,'.');
    progress = n/reps;
    waitbar(progress,pb,sprintf('%.0f%%',progress*100));

    % Initial number of reactant molecules; Assume B is limiting for 2nd order rx.
    Ao = 100;                       Bo = floor(0.85*Ao);      
    Co = floor(1.2*Ao);             Do = floor(0.75*Ao);   
    % Arrays for tracking reactant agent status
    A = ones(1,Ao);                 B = ones(1,Bo);          
    C = ones(1,Co);                 D = ones(1,Do);          
    % Initialize time-dependent sum of molecule numbers
    At = zeros(1,steps);            Bt = zeros(1,steps);     
    Ct = zeros(1,steps);            Dt = zeros(1,steps);
    At(1) = sum(A);                 Bt(1) = sum(B);          
    Ct(1) = sum(C);                 Dt(1) = sum(D);        
    Xt = zeros(1,steps);            Yt = zeros(1,steps);
    Xt(1) = 0;                      Yt(1) = 0;

    time = zeros(1,steps);

    for t=2:steps                   %  while Bt(t-1) > 0 && At(t-1) > 0
        
        deltaX = 0;             deltaY = 0;  

    %     dt = 10 / At(t-1);                        % Variable time step incr. 1
    %     dt = 1/(k1*At(t-1) + k2*At(t-1)*Bt(t-1)); % Variable time step incr. 2
    %     dt = exprnd(1/(k1*At(t-1) + k2*At(t-1)*Bt(t-1)));   % Exponentially-distributed time step increment
        
%         P1 = k1 * Bt(t-1) * dt;                   % P_dif of 2nd order rx A+B-->X
        P1 = 1 - exp(-k1 * Bt(t-1) * dt);         % P_ber of 2nd order rx A+B-->X
%         P2 = k2 * Ct(t-1) * Dt(t-1) * dt;         % P_dif-A of 3rd order rx A+C+D-->Y
        P2 = 1 - exp(-k2 * Ct(t-1) * Dt(t-1) * dt);   % P_ber-A of 3rd order rx A+C+D-->Y
        
        tempA = find(A==1);    
        for i = 1:size(tempA,2)
            r = rand;
            if r < P1                           % 2nd order rx
                tempB = find(B==1);
                b = ceil(rand*size(tempB,2));
                if b==0,            continue;           end
                A(tempA(i)) = 0;  
                B(tempB(b)) = 0;
                deltaX = deltaX + 1;
            elseif r >= P1 && r < P1+P2         % 3rd order rx
                tempB = find(B==1);      tempC = find(C==1);   tempD = find(D==1);
                b = ceil(rand*size(tempB,2));
                c = ceil(rand*size(tempC,2));
                d = ceil(rand*size(tempD,2));
                if b==0 || c==0 || d==0,         continue;           end
                A(tempA(i)) = 0;
                C(tempC(c)) = 0;
                D(tempD(d)) = 0;
                deltaY = deltaY + 1;
            end 
        end 
        At(t) = sum(A);                 Bt(t) = sum(B);        
        Ct(t) = sum(C);                 Dt(t) = sum(D);
        Xt(t) = Xt(t-1) + deltaX;       Yt(t) = Yt(t-1) + deltaY;
        time(t) = time(t-1) + dt;
    end
    Sa(n,:) = At;             Sb(n,:) = Bt;          
    Sc(n,:) = Ct;             Sd(n,:) = Dt;
    Sx(n,:) = Xt;             Sy(n,:) = Yt;

end         % end 'for n' loop

avgA = mean(Sa);                     sdevA = std(Sa);
avgB = mean(Sb);                     sdevB = std(Sb);
avgC = mean(Sc);                     sdevC = std(Sc);
avgD = mean(Sd);                     sdevD = std(Sd);
avgX = mean(Sx);                     sdevX = std(Sx);
avgY = mean(Sy);                     sdevY = std(Sy);

%% Solve ODE for mixed order kinetics
% Same time sampling as in simulation
[t_sol,y_sol] = ode45(@o23_mixed2_dif,time,[At(1);Bt(1);Ct(1);Dt(1);Xt(1);Yt(1)]);
% Or, uncomment the following to try a different time sampling
% [t_sol,y_sol] = ode45(@o23_mixed2_dif,0:t_max/500:t_max,[At(1);Bt(1);Ct(1);Dt(1);Xt(1);Yt(1)]);

%% Find Coefficient of Determination, R^2
% This works for any curve fit, linear or nonlinear.
SST(1) = sum((avgA' - mean(avgA)).^2);     % Total sum of squares for simulation data (for A)
SSR(1) = sum((avgA' - y_sol(:,1)).^2);     % sum of square residuals (sim data vs DE predictions)

SST(2) = sum((avgB' - mean(avgB)).^2);     % Total sum of squares for simulation data (for B)
SSR(2) = sum((avgB' - y_sol(:,2)).^2);     % sum of square residuals (sim data vs DE predictions)

SST(3) = sum((avgC' - mean(avgC)).^2);     % Total sum of squares for simulation data (for C)
SSR(3) = sum((avgC' - y_sol(:,3)).^2);     % sum of square residuals (sim data vs DE predictions)

SST(4) = sum((avgD' - mean(avgD)).^2);     % Total sum of squares for simulation data (for D)
SSR(4) = sum((avgD' - y_sol(:,4)).^2);     % sum of square residuals (sim data vs DE predictions)

SST(5) = sum((avgX' - mean(avgX)).^2);     % Total sum of squares for simulation data (for X)
SSR(5) = sum((avgX' - y_sol(:,5)).^2);     % sum of square residuals (sim data vs DE predictions)

SST(6) = sum((avgY' - mean(avgY)).^2);     % Total sum of squares for simulation data (for Y)
SSR(6) = sum((avgY' - y_sol(:,6)).^2);     % sum of square residuals (sim data vs DE predictions)


Rsq = 1 - SSR./SST                         % Definition of R^2
R = sqrt(Rsq);                             % Correlation Coefficient R

%% Plot Time Courses
figure1 = figure('Name','Mixed o2 Rx Time course','NumberTitle','off');
set(figure1,'Position',[1 1 500 450]);
p1 = plot(time,avgA,'b','MarkerSize',3,'DisplayName','<N_A(t)>_{sim}');                 hold on;
p1_dev1 = plot(time,avgA+sdevA,'LineStyle','--','Color',[0.8 0.8 0.8]);
p1_dev0 = plot(time,avgA-sdevA,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p1s = plot(time,Sa(end,:),'m','MarkerSize',3,'DisplayName','Sample Trajectory');

p2 = plot(time,avgB,'r','MarkerSize',3,'DisplayName','<N_B(t)>_{sim}');                 
% p2_dev1 = plot(time,avgB+sdevB,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p2_dev0 = plot(time,avgB-sdevB,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p2s = plot(time,Sb(end,:),'m','MarkerSize',3,'DisplayName','Sample Trajectory');

p3 = plot(time,avgC,'c','MarkerSize',3,'DisplayName','<N_C(t)>_{sim}');                 
% p3_dev1 = plot(time,avgC+sdevC,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p3_dev0 = plot(time,avgC-sdevC,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p3s = plot(time,Sc(end,:),'m','MarkerSize',3,'DisplayName','Sample Trajectory');

p4 = plot(time,avgD,'Color',[0.9 0.9 0],'MarkerSize',3,'DisplayName','<N_D(t)>_{sim}');                 
% p4_dev1 = plot(time,avgD+sdevD,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p4_dev0 = plot(time,avgD-sdevD,'LineStyle','--','Color',[0.8 0.8 0.8]);
p4s = plot(time,Sd(end-3,:),'m','MarkerSize',3,'DisplayName','Sample Trajectory');

p5 = plot(time,avgX,'k','MarkerSize',3,'DisplayName','<N_X(t)>_{sim}');                 
% p5_dev1 = plot(time,avgX+sdevX,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p5_dev0 = plot(time,avgX-sdevX,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p5s = plot(time,Sx(end-16,:),'m','MarkerSize',3,'DisplayName','Sample Trajectory');

p6 = plot(time,avgY,'g','MarkerSize',3,'DisplayName','<N_Y(t)>_{sim}');                 
p6_dev1 = plot(time,avgY+sdevY,'LineStyle','--','Color',[0.8 0.8 0.8]);
p6_dev0 = plot(time,avgY-sdevY,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p6s = plot(time,Sy(end,:),'m','MarkerSize',3,'DisplayName','Sample Trajectory');

p1d = plot(t_sol,y_sol(:,1),':b','DisplayName','DE');        
p2d = plot(t_sol,y_sol(:,2),':r');
p3d = plot(t_sol,y_sol(:,3),':c');
p4d = plot(t_sol,y_sol(:,4),':','Color',[0.9 0.9 0]);
p5d = plot(t_sol,y_sol(:,5),':k');
p5d = plot(t_sol,y_sol(:,6),':g');

xlabel('t (sec)');                 ylabel('N(t)');                  hold off;   
axis tight;                         % axis([0 t_max 0 agents]);                                           
title({['A + B        \rightarrow X  (k_1 = ' num2str(k1) '0 sec^{-1})'],...
    ['A + C + D \rightarrow Y  (k_2 = ' num2str(k2) ' sec^{-1})']}...
    ,'FontName','Times New Roman','FontSize',11)
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
leg = legend([p1 p2 p3 p4 p5 p6]);
% set(leg,'OuterPosition',[0.673 0.374 0.182 0.123]);
set(leg,'Location','NorthEast');
set(leg,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95]);
% set(p??,'Visible','off');

%% Finish
clear temp x deltaC deltaD r i;
toc

% Result: It works!
% - Works for fixed time step increments.
% - Note that calculated probabilities work best for small dt. For variable
% time step increments some dt's will be large (towards the end of the
% simulation). It may be best to use fixed time step increments while
% calculating such probabilities, despite its computational cost.