% 2A + B --> D
% Algorithm is run wrt A.
% Simulating 3rd order kinetics using my agent-based algorithm.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear; clc;                         tic

global k;       % Termolecular microscopic kinetic constant [units: 1/sec];
k = 0.01;

agents = 15;
reps = 100;

t_max = 10;                             % in seconds
dt = 1/100;                             % Fixed time step increment
steps = t_max / dt;

Sa = zeros(reps,steps);      Sb = zeros(reps,steps);      Sd = zeros(reps,steps);

% Initial number of reactant molecules;
Ao = agents;      Bo = floor(0.4*agents);          Dmax = Bo;

for n=1:reps
    % Arrays for tracking reactant and product agent status
    A = ones(1,Ao);        B = ones(1,Bo);      D = zeros(1,Dmax); 
    % Initialize time-dependent sum of molecule numbers
    At(1) = sum(A);        Bt(1) = sum(B);      Dt(1) = 0;
    
    time = zeros(1,steps);
    
    for t=2:steps                     
    %     dt = 10  / ( (At(t-1)-1) * Bt(t-1) );       % Variable time step increment

    %     P = k * (At(t-1)-1) * Bt(t-1) * dt;
        P = 1 - exp(-k * (At(t-1)-1) * Bt(t-1) * dt);

        for i = 1:size(A,2)
            if A(i) == 1 && rand < P
                A(i) = 0;                       % First Agent A dies
                temp = find(A==1);
                w = ceil(rand*size(temp,2));    % Pick another A agent randomly
                A(temp(w)) = 0;                 % Second agent A dies
                temp2 = find(B==1);      
                x = ceil(rand*size(temp2,2));   % Pick B agent randomly
                B(temp2(x)) = 0;                % agent B dies
                
                D(temp2(x)) = 1;                % Agent D is born!
            end
        end 
        time(t) = time(t-1) + dt;
        At(t) = sum(A);         Bt(t) = sum(B);         Dt(t) = sum(D);
    end
    Sa(n,:) = At;           Sb(n,:) = Bt;                Sd(n,:) = Dt;

end
avgA = mean(Sa);                     sdevA = std(Sa);
avgB = mean(Sb);                     sdevB = std(Sb);
avgD = mean(Sd);                     sdevD = std(Sd);

%% Solve ODE for 3rd order kinetics
[t_sol,y_sol] = ode45(@o3_het_dif,[0:t_max/200:t_max],[At(1) ; Bt(1) ; Dt(1)]);

%% Plot
figure1 = figure('Name','3rd Order Rx Time course','NumberTitle','off');
set(figure1,'Position',[1 1 500 406]);
p1 = plot(time,avgA,'b','MarkerSize',3,'DisplayName','A');                 hold on;
p1_dev1 = plot(time,avgA+sdevA,'LineStyle','--','Color',[0.8 0.8 0.8]);
p1_dev0 = plot(time,avgA-sdevA,'LineStyle','--','Color',[0.8 0.8 0.8]);
p2 = plot(time,avgB,'r','MarkerSize',3,'DisplayName','B');                 
p2_dev1 = plot(time,avgB+sdevB,'LineStyle','--','Color',[0.8 0.8 0.8]);
p2_dev0 = plot(time,avgB-sdevB,'LineStyle','--','Color',[0.8 0.8 0.8]);
p3 = plot(time,avgD,'g','MarkerSize',3,'DisplayName','D');                 
p3_dev1 = plot(time,avgD+sdevD,'LineStyle','--','Color',[0.8 0.8 0.8]);
p3_dev0 = plot(time,avgD-sdevD,'LineStyle','--','Color',[0.8 0.8 0.8]);

p1d = plot(t_sol,y_sol(:,1),':b','DisplayName','DE_{ab}');        
p2d = plot(t_sol,y_sol(:,2),':r');
p3d = plot(t_sol,y_sol(:,3),':g');

xlabel('t (sec)');                 ylabel('N(t)');           
axis([0 t_max 0 agents]);                                           hold off;
title(['N_{A,i} = ' num2str(Ao) ' ; N_{B,i} = ' num2str(Bo)],...
    'FontName','Times New Roman','FontSize',12)
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
leg = legend([p1 p2 p3]);
% set(leg,'OuterPosition',[0.673 0.374 0.182 0.123]);
set(leg,'Location','NorthEast');
set(leg,'FontName','Times New Roman','EdgeColor',[0.95 0.95 0.95]);
% set(p??,'Visible','off');

%% Finish
clear temp temp2 x i n p*;
% save(['o3_A=' num2str(agents) '.mat']);

toc

% Result:
% - It works for both fixed and variable time step increment algorithms. 
% - It's always best to run the algorithm wrt the limiting reactant!
%   (If done wrt xs reactant, complications arise during the algorithm's 
%   implementation when the pool of limiting reactant has been exhausted.
%   Can be dealt with, but at the cost of programmatic simplicity.)
