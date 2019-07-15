%  A + B + C --> D
% Simulating 3rd order kinetics using my agent-based algorithm.
% Assume volume is --- microliters.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear; clc;                 tic

global k;       % Termolecular microscopic kinetic constant [units: 1/sec];
k = 0.001;

%% Initial number of reactant molecules; assume B is limiting.
agents = 10;
reps = 500;

t_max = 50;                             % in seconds
dt = 1/100;                             % Fixed time step increment
steps = t_max / dt;

Sa = zeros(reps,steps);      Sb = zeros(reps,steps);        Sc = zeros(reps,steps);      
Sd = zeros(reps,steps);

% Initial number of reactant molecules;
Ao = agents;      Bo = 0.7*agents;       Co = 0.90*agents;       Dmax = Bo;

for n=1:reps

    % Arrays for tracking reactant and product agent status
    A = ones(1,Ao);        B = ones(1,Bo);      C = ones(1,Co);          
    D = zeros(1,Dmax); 
    % Initialize time-dependent sum of molecule numbers
    At(1) = sum(A);        Bt(1) = sum(B);      Ct(1) = sum(C);     
    Dt(1) = 0;

    time = zeros(1,steps);

    for t=2:steps                     
    %     dt = 10  / ( At(t-1)*Ct(t-1) );       % Variable time step increment

        P = k * At(t-1) * Ct(t-1)  * dt;
%         P = 1 - exp(-k * At(t-1) * Ct(t-1)  * dt);
        
        temp = find(B==1);
        for i = 1:size(temp,2)
            if B(temp(i)) == 1 && rand < P
                B(temp(i)) = 0;           % Agent B dies
                A(temp(i)) = 0;           % Agent A dies
                C(temp(i)) = 0;           % Agent C dies
                D(temp(i)) = 1;           % Agent D is born!
            end
        end 
        time(t) = time(t-1) + dt;
        At(t) = sum(A);         Bt(t) = sum(B);         Ct(t) = sum(C);
        Dt(t) = sum(D);
    end
    Sa(n,:) = At;           Sb(n,:) = Bt;               Sc(n,:) = Ct;       
    Sd(n,:) = Dt;

end
avgA = mean(Sa);                     sdevA = std(Sa);
avgB = mean(Sb);                     sdevB = std(Sb);
avgC = mean(Sc);                     sdevC = std(Sc);
avgD = mean(Sd);                     sdevD = std(Sd);

%% Solve ODE for 3rd order kinetics - Plot Time Courses
[t_sol, y_sol] = ode45(@o3_het_dif,[0:t_max/200:t_max],[At(1) ; Bt(1) ; Ct(1) ; Dt(1)]);

%% Plot
figure1 = figure('Name','3rd Order Rx Time course','NumberTitle','off');
set(figure1,'Position',[1 1 500 406]);
p1 = plot(time,avgA,'b','MarkerSize',3,'DisplayName','A');                 hold on;
p1_dev1 = plot(time,avgA+sdevA,'LineStyle','--','Color',[0.8 0.8 0.8]);
p1_dev0 = plot(time,avgA-sdevA,'LineStyle','--','Color',[0.8 0.8 0.8]);
p2 = plot(time,avgB,'r','MarkerSize',3,'DisplayName','B');                 
p2_dev1 = plot(time,avgB+sdevB,'LineStyle','--','Color',[0.8 0.8 0.8]);
p2_dev0 = plot(time,avgB-sdevB,'LineStyle','--','Color',[0.8 0.8 0.8]);
p3 = plot(time,avgC,'k','MarkerSize',3,'DisplayName','C');                 
p3_dev1 = plot(time,avgC+sdevC,'LineStyle','--','Color',[0.8 0.8 0.8]);
p3_dev0 = plot(time,avgC-sdevC,'LineStyle','--','Color',[0.8 0.8 0.8]);
p4 = plot(time,avgD,'g','MarkerSize',3,'DisplayName','D');                 
p4_dev1 = plot(time,avgD+sdevD,'LineStyle','--','Color',[0.8 0.8 0.8]);
p4_dev0 = plot(time,avgD-sdevD,'LineStyle','--','Color',[0.8 0.8 0.8]);

p1d = plot(t_sol,y_sol(:,1),':b','DisplayName','DE_{ab}');        
p2d = plot(t_sol,y_sol(:,2),':r');
p3d = plot(t_sol,y_sol(:,3),':k');
p4d = plot(t_sol,y_sol(:,4),':g');

xlabel('t (sec)');                 ylabel('N(t)');           
axis([0 t_max 0 agents]);                                           hold off;
title(['N_{A,i} = ' num2str(Ao) ' ; N_{B,i} = ' num2str(Bo) ' ; N_{C,i} = ' num2str(Co)],...
    'FontName','Times New Roman','FontSize',12)
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
leg = legend([p1 p2 p3 p4]);
% set(leg,'OuterPosition',[0.673 0.374 0.182 0.123]);
set(leg,'Location','NorthEast');
set(leg,'FontName','Times New Roman','EdgeColor',[0.95 0.95 0.95]);
% set(p??,'Visible','off');

%% Finish
clear temp x i p*;
toc

% Result: 
% - It works for both fixed and variable time step increment algorithms. 
