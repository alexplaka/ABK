%  3A --> D
% Simulating 3rd order kinetics using my agent-based algorithm.
% The integrated rate law for this rx is:
% A/Ao = sqrt( 1 / (2 * 3 * k * Ao^2 * t + 1))
% This is used to calculate the probability of rx for time duration dt.
% Notice that the stoichiometric coeeficient, 3, is omitted from the
% the expression for P, since P is evaluated for each molecule
% separately from the others, and the correct number of molecules is
% marked as reacted when rand < P.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;  clc;                                 tic
rng(0);

global k dif_canon;                       

agents = 12;

% Avonum = 6.02e+23;                    % Avogadro's number
% V = 10^-21;                           % Volume in Liters
% k = km /(Avonum*V)^2;                 % Bimolecular molar kinetic constant [units: 1 /(M sec)];
k = 0.001;                               

t_max = 50;                             % in seconds
dt = 1/100;                             % Fixed time step increment
steps = t_max / dt;
% dt = 1 / (20*k*Ao^2);                   % Fixed time step increment
% The above equation for dt ensures that the initial probability is <1
% The factor of 20 is arbitrary and assures that P(1)<0.05
time = zeros(1,steps);

reps = 500;
Sa = zeros(reps,steps);                 Sd = zeros(reps,steps);
P = zeros(reps,steps);
%% ABK algorithm
for n=1:reps
    
    % Initial number of reactant molecules; assume B is limiting.
    Ao = agents;                    Dmax = Ao;
    % Arrays for tracking reactant and product agent status
    A = ones(1,Ao);                 D = zeros(1,Dmax); 
    % Initialize time-dependent sum of molecule numbers
    At = zeros(1,steps);            Dt = zeros(1,steps);
    At(1) = sum(A);                 Dt(1) = sum(D);

    time(1) = 0;                    

    for t=2:steps                     
%           dt = 10  / ( (At(t-1)-1) * (At(t-1)-2) );       % Variable time step increment

%         P(n,t-1) = 1 - sqrt(1 / (2 * k * At(t-1)^2 * dt + 1));                  % P_can
%         P(n,t-1) =   DOES NOT EXIST!!!                                          % P_int
        P(n,t-1) = 1 - exp(-k * (At(t-1)-1) * (At(t-1)-2) * dt);                % P_ber
%         P(n,t-1) = k * (At(t-1)-1) * (At(t-1)-2) * dt;                          % P_dif

        temp = find(A==1);
        if size(temp,2) <= 2          
            At(t:end) = At(t-1);
            break;          
        end
                
        for i = 1:size(temp,2)
            if A(temp(i))==1 && rand < P(n,t-1)
            % First IF condition ensures that an A agent has not already "died" in this time step
                A(temp(i)) = 0;           % First agent A dies
                temp2 = find(A==1);      
                x = ceil(rand(1,2)*size(temp2,2));  % Pick 2 A agents randomly
                A(temp2(x(1))) = 0;        % Second agent A dies
                A(temp2(x(2))) = 0;        % Third agent A dies
                D(temp(i)) = 1;            % Agent D is born!
            end
        end 

        At(t) = sum(A);                 Dt(t) = sum(D);
        time(t) = time(t-1) + dt;
    end
    Sa(n,:) = At;                       Sd(n,:) = Dt;

end         % end 'for n' loop
avg = mean(Sa);                     sdev = std(Sa);

%% Do the same using Gillespie's algorithm
% gA(1) =  agents;            gC(1) = 0;
% for n=1:reps
%     gtime(1) = 0;
%     t = 2;
%     while gA(t-1) > 1                                
%         a(t) = k * gA(t-1) * (gA(t-1)-1) * (gA(t-1)-2);                   % Propensity function
%         gdt = -log(rand) / a(t);
%         gtime(t) = gtime(t-1) + gdt;
%         gA(t) = gA(t-1) - 3;
%         gD(t) = gD(t-1) + 1;    
%         t = t + 1;
%     end
%     gtime_all(n,:) = gtime;
%     gA_all(n,:) = gA;                   gD_all(n,:) = gD;
% end
% avg_gtime = mean(gtime_all);
%% Solve ODE for 2nd order kinetics
dif_canon = 0;                  % use agent-based form of DE
[t_sol0, y_sol0] = ode45(@o3_homo_dif,0:t_max/200:t_max,[At(1) ; Dt(1)]);
dif_canon = 1;                  % use canonical form of DE
[t_sol1, y_sol1] = ode45(@o3_homo_dif,0:t_max/200:t_max,[At(1) ; Dt(1)]);
%% Plot
figure1 = figure('Name','3rd Order Rx Time course','NumberTitle','off');
set(figure1,'Position',[1 1 500 406]);
p1 = plot(time,avg,'b','MarkerSize',3,'DisplayName','ABK');                         hold on;
p1_dev1 = plot(time,avg+sdev,'LineStyle','--','Color',[0.8 0.8 0.8]);
p1_dev0 = plot(time,avg-sdev,'LineStyle','--','Color',[0.8 0.8 0.8]);
xlabel('t (sec)');                 ylabel('N_A(t)');           
p2 = plot(t_sol0,y_sol0(:,1),':g','DisplayName','DE_{ab}');        % plot(t_sol0,y_sol0(:,2),'r');
p3 = plot(t_sol1,y_sol1(:,1),':r','DisplayName','DE_{can}');
% p4 = plot(avg_gtime',gA_all(1,:)','LineStyle',':','Color','k','DisplayName','Gil'); % Gillespie trajectories
axis([0 t_max 0 agents]);                                           hold off;
title(['N_{A,i} = ' num2str(agents)],'FontName','Times New Roman','FontSize',12)
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
leg = legend([p1 p2 p3]);
% set(leg,'OuterPosition',[0.673 0.374 0.182 0.123]);
set(leg,'Location','NorthEast');
set(leg,'FontName','Times New Roman','EdgeColor',[0.95 0.95 0.95]);
% set(p2,'Visible','off');                        set(p3,'Visible','off');

%% Finish
toc
clear temp temp2 x i tmin tmax idxmin idxmax p*;
% save(['o2_A=' num2str(agents) '.mat']);
