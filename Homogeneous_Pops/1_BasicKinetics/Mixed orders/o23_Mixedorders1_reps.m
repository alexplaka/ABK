%  A + B     --> X            Rate constant k1  (2nd order)
%  A + B + C --> Y            Rate constant k2  (3rd order)
% Simulating concurrent 2nd and 3rd order kinetics using my agent-based
% algorithm.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear; clc;  tic
rng(0);
pb = waitbar(0,'0');

global k1 k2;
k1 = 0.05;              % k1: 2nd order rate constant [units: 1/sec]
k2 = 0.001;             % k2: 3rd order *microscopic* rate constant [units: 1/sec]

t_max = 5;                             % in seconds
dt = 1/500;                            % Fixed time step increment
steps = t_max / dt;

reps = 500;
Sa = zeros(reps,steps);     Sb = zeros(reps,steps);     Sc = zeros(reps,steps);
Sx = zeros(reps,steps);     Sy = zeros(reps,steps);

P1a = zeros(reps,steps);     P2a = zeros(reps,steps);

for n=1:reps
%     fprintf(1,'.');
    progress = n/reps;
    waitbar(progress,pb,sprintf('%.0f%%',progress*100));

    % Initial number of reactant molecules; Assume B is limiting for 2nd order rx.
    Ao = 10;              Bo = floor(0.7*Ao);      Co = floor(1.2*Ao);   
    % Arrays for tracking reactant agent status
    A = ones(1,Ao);        B = ones(1,Bo);          C = ones(1,Co);          
    % Initialize time-dependent sum of molecule numbers
    At = zeros(1,steps);   Bt = zeros(1,steps);     Ct = zeros(1,steps);
    At(1) = sum(A);        Bt(1) = sum(B);          Ct(1) = sum(C);        
    Xt = zeros(1,steps);            Yt = zeros(1,steps);
    Xt(1) = 0;                      Yt(1) = 0;

    time = zeros(1,steps);

    for t=2:steps                   %  while Bt(t-1) > 0 && At(t-1) > 0
        
        deltaX = 0;             deltaY = 0;  

    %     dt = 10 / At(t-1);                        % Variable time step incr. 1
    %     dt = 1/(k1*At(t-1) + k2*At(t-1)*Bt(t-1)); % Variable time step incr. 2
    %     dt = exprnd(1/(k1*At(t-1) + k2*At(t-1)*Bt(t-1)));   % Exponentially-distributed time step increment
        
        P1a(n,t-1) = k1 * Bt(t-1) * dt;                   % P_dif of 2nd order rx A+B-->X
%         P1a(n,t-1) = 1 - exp(-k1 * Bt(t-1) * dt);         % P_ber of 2nd order rx A+B-->X
        P2a(n,t-1) = k2 * Bt(t-1) * Ct(t-1) * dt;         % P_dif-A of 3rd order rx A+B+C-->Y
%         P2a(n,t-1) = 1 - exp(-k2 * Bt(t-1) * Ct(t-1) * dt);   % P_ber-A of 3rd order rx A+B+C-->Y
        
%         P1b = k1 * At(t-1) * dt;                   % P_dif of 2nd order rx A+B-->X
%         P1b = 1 - exp(-k1 * At(t-1) * dt);         % P_ber of 2nd order rx A+B-->X
%         P2b = k2 * At(t-1) * Ct(t-1) * dt;         % P_dif-B of 3rd order rx A+B+C-->Y
%         P2b = 1 - exp(-k2 * At(t-1) * Ct(t-1) * dt);   % P_ber-B of 3rd order rx A+B+C-->Y

        % *** Model two processes at the same time wrt A! ***
        tempA = find(A==1);     
        
        for i = 1:size(tempA,2)
            r = rand; 
            if r < P1a(n,t-1)                          % if 2nd order rx wrt A
                tempB = find(B==1);
                b = ceil(rand*size(tempB,2));
                if b==0,            continue;           end
                A(tempA(i)) = 0;  
                B(tempB(b)) = 0;         % tempB(b)=[];   
                deltaX = deltaX + 1;
            elseif r >= P1a(n,t-1) && r < P1a(n,t-1)+P2a(n,t-1)      % 3rd order rx wrt A
                tempB = find(B==1);   tempC = find(C==1);
                b = ceil(rand*size(tempB,2));
                c = ceil(rand*size(tempC,2));
                if b==0 || c==0,         continue;           end
                A(tempA(i)) = 0;             
                B(tempB(b)) = 0;         % tempB(b)=[];    
                C(tempC(c)) = 0;         % tempC(c)=[];
                deltaY = deltaY + 1;
            end
        end
                      
        At(t) = sum(A);          Bt(t) = sum(B);        Ct(t) = sum(C);
        Xt(t) = Xt(t-1) + deltaX;       Yt(t) = Yt(t-1) + deltaY;
        time(t) = time(t-1) + dt;

    end     % "for t"
                
    Sa(n,:) = At;             Sb(n,:) = Bt;          Sc(n,:) = Ct;
    Sx(n,:) = Xt;             Sy(n,:) = Yt;

end         % end 'for n' loop

avgA = mean(Sa);                     sdevA = std(Sa);
avgB = mean(Sb);                     sdevB = std(Sb);
avgC = mean(Sc);                     sdevC = std(Sc);
avgX = mean(Sx);                     sdevX = std(Sx);
avgY = mean(Sy);                     sdevY = std(Sy);

%% Solve ODE for mixed order kinetics
% Same time sampling as in simulation
[t_sol,y_sol] = ode45(@o23_mixed1_dif,time,[At(1);Bt(1);Ct(1);Xt(1);Yt(1)]);
% Or, uncomment the following to try a different time sampling
% [t_sol,y_sol] = ode45(@o23_mixed1_dif,0:t_max/500:t_max,[At(1);Bt(1);Ct(1);Xt(1);Yt(1)]);

%% Find Coefficient of Determination, R^2
% This works for any curve fit, linear or nonlinear.
SST(1) = sum((avgA' - mean(avgA)).^2);     % Total sum of squares for simulation data (for A)
SSR(1) = sum((avgA' - y_sol(:,1)).^2);     % sum of square residuals (sim data vs DE predictions)

SST(2) = sum((avgB' - mean(avgB)).^2);     % Total sum of squares for simulation data (for B)
SSR(2) = sum((avgB' - y_sol(:,2)).^2);     % sum of square residuals (sim data vs DE predictions)

SST(3) = sum((avgC' - mean(avgC)).^2);     % Total sum of squares for simulation data (for C)
SSR(3) = sum((avgC' - y_sol(:,3)).^2);     % sum of square residuals (sim data vs DE predictions)

Rsq = 1 - SSR./SST                         % Definition of R^2
R = sqrt(Rsq);                             % Correlation Coefficient R

%% Plot Time Courses
figure1 = figure('Name','Mixed o2 Rx Time course','NumberTitle','off');
set(figure1,'Position',[1 1 500 450]);
p1 = plot(time,avgA,'b','MarkerSize',3,'DisplayName','<N_A(t)>_{sim}');                 hold on;
% p1_dev1 = plot(time,avgA+sdevA,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p1_dev0 = plot(time,avgA-sdevA,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p1s = plot(time,Sa(end,:),'m','MarkerSize',3,'DisplayName','Sample Trajectory');

p2 = plot(time,avgB,'r','MarkerSize',3,'DisplayName','<N_B(t)>_{sim}');                 
% p2_dev1 = plot(time,avgB+sdevB,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p2_dev0 = plot(time,avgB-sdevB,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p2s = plot(time,Sb(end,:),'m','MarkerSize',3,'DisplayName','Sample Trajectory');

p3 = plot(time,avgC,'c','MarkerSize',3,'DisplayName','<N_C(t)>_{sim}');                 
% p3_dev1 = plot(time,avgC+sdevC,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p3_dev0 = plot(time,avgC-sdevC,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p3s = plot(time,Sc(end,:),'m','MarkerSize',3,'DisplayName','Sample Trajectory');

p4 = plot(time,avgX,'k','MarkerSize',3,'DisplayName','<N_X(t)>_{sim}');                 
% p4_dev1 = plot(time,avgX+sdevX,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p4_dev0 = plot(time,avgX-sdevX,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p4s = plot(time,Sx(end-16,:),'m','MarkerSize',3,'DisplayName','Sample Trajectory');

p5 = plot(time,avgY,'g','MarkerSize',3,'DisplayName','<N_Y(t)>_{sim}');                 
% p5_dev1 = plot(time,avgY+sdevY,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p5_dev0 = plot(time,avgY-sdevY,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p5s = plot(time,Sy(end,:),'m','MarkerSize',3,'DisplayName','Sample Trajectory');

p1d = plot(t_sol,y_sol(:,1),':b','DisplayName','DE');        
p2d = plot(t_sol,y_sol(:,2),':r');
p3d = plot(t_sol,y_sol(:,3),':c');
p4d = plot(t_sol,y_sol(:,4),':k');
p5d = plot(t_sol,y_sol(:,5),':g');

xlabel('t (sec)');                 ylabel('N(t)');                  hold off;   
axis tight;                         % axis([0 t_max 0 agents]);                                           
title({['A + B     \rightarrow X  (k_1 = ' num2str(k1) ' sec^{-1})'],...
    ['A + B + C \rightarrow Y  (k_2 = ' num2str(k2) ' sec^{-1})']}...
    ,'FontName','Times New Roman','FontSize',11)
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
leg = legend([p1 p2 p3 p4 p5]);
% set(leg,'OuterPosition',[0.673 0.374 0.182 0.123]);
set(leg,'Location','NorthEast');
set(leg,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95]);
% set(p??,'Visible','off');

%% Finish
close(pb);
clear progress pb temp x deltaC deltaD r i p* fig* leg;
toc
