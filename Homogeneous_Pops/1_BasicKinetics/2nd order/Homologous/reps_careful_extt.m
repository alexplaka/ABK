%  A + A --> C
% Simulating 2nd order kinetics using my agent-based algorithm.
% The integrated rate law for this rx is:
% A/Ao = 1 / (2 * k * Ao * t + 1)
% This is used to calculate the probability of rx for time duration dt.
% Notice that the stoichiometric coeeficient, 2, is omitted from the
% the expression for P, since P is evaluated for *each* molecule
% separately from the others, and the correct number of molecules is
% marked as reacted when rand < P.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

clear;                                  tic
rng(0);

global k dif_canon;                       

agents = 10;

% Initial number of agents
Ao = agents;                    Co = 0;                 Cmax = Ao;

% Avonum = 6.02e+23;                    % Avogadro's number
% V = 10^-21;                           % Volume in Liters
% km = 0.1;                     % Bimolecular molar kinetic constant [units: 1 /(M sec)];
k = 0.01;                               % k = km /(Avonum*V);

reps = 2000;

t_max = 30;                           % in seconds
dt = 1/100;                              % Fixed time step increment
steps = t_max / dt;
time = zeros(1,steps);

Sa = zeros(reps,steps);                 Sc = zeros(reps,steps);
P = zeros(reps,steps);
extt = zeros(reps, Ao);             % *** To monitor extinction time of A agents **

%% ABK algorithm
for n=1:reps
    
    fprintf('.');
    % Arrays for tracking reactant and product agent status
    A = ones(1,Ao);                 C = zeros(1,Cmax); 
    % Initialize time-dependent sum of molecule numbers
    At = zeros(1,steps);            Ct = zeros(1,steps);
    At(1) = sum(A);                 Ct(1) = sum(C);

    time(1) = 0;

    for t=2:steps                     
    %     dt = 10  / At(t-1);             % Variable time step increment
    %     dt = 1/( k*At(t-1) );   % Variable time step increment (ver. 2)
    %     dt = exprnd(1/( k*At(t-1) ) )   % Exponentially-distributed time step increment
%       if At(t-1) > 1
%         P(n,t) = 1 - 1 / (1 + k * (At(t-1)-1) * dt);             % P version 1  (P_can)
%         P(n,t) = 1 - 1 /( At(t-1)*(1-exp(-k*dt)) + exp(-k*dt) ); % P version 2  (P_int)
        P(n,t) = 1 - exp(-k * (At(t-1)-1) * dt);                 % P version 3  (P_ber)
%         P(n,t) = k * (At(t-1)-1) * dt;                           % P version 4 (P_dif)

        tempAa = find(A==1);                % Find alive A agents
        for i = 1:size(tempAa,2)
            if tempAa(i) ~= 0  && rand < P(n,t) 
                A(tempAa(i)) = 0;                   % First A agent dies
                
                x = ceil(rand*size(tempAa,2));
                while tempAa(x) == 0 || x == i       
                    x = ceil(rand*size(tempAa,2));
                end
                A(tempAa(x)) = 0;             % Second A agent dies (chosen randomly)
                
                C(i) = 1;   
                
                extt(n,tempAa(i)) = time(t-1) + dt;
                extt(n,tempAa(x)) = time(t-1) + dt;
                
                tempAa(i) = 0;       % Make sure this A agent is not reselected in this time step
                tempAa(x) = 0;       % Make sure this A agent is not reselected in this time step
            end
        end
        
        At(t) = sum(A);                 Ct(t) = sum(C);
        time(t) = time(t-1) + dt;
    end
    
    Sa(n,:) = At;                       Sc(n,:) = Ct;

end         % end 'for n' loop

avg = mean(Sa);                     sdev = std(Sa);
% avgC = mean(Sc);                     sdevC = std(Sc);

%% Remove zeros from extt (zeros represent agents that didn't "die" within totalTime)
f = size(extt,2);                      
for z=1:f
    extt_zeros = find(extt(:,z)==0);
    for a = 1:size(extt_zeros,1)
        extt(extt_zeros(a),z) = NaN;
    end
end

extt = reshape(extt,[],1);            % Reshape to column vector

%% Do the same using Gillespie's algorithm
% gA(1) =  agents;            gC(1) = 0;
% for n=1:reps
%     gtime(1) = 0;
%     t = 2;
%     while gA(t-1) > 1                                
%         a(t) = k * gA(t-1) * (gA(t-1)-1);                         % Propensity function
%         gdt = -log(rand) / a(t);
%         gtime(t) = gtime(t-1) + gdt;
%         gA(t) = gA(t-1) - 2;
%         gC(t) = gC(t-1) + 1;    
%         t = t + 1;
%     end
%     gtime_all(n,:) = gtime;
%     gA_all(n,:) = gA;                   gC_all(n,:) = gC;
% end
% avg_gtime = mean(gtime_all);
%% Solve ODE for 2nd order kinetics
dif_canon = 0;                  % use agent-based form of DE
[t_sol0, y_sol0] = ode45(@o2_homo_dif,0:t_max/200:t_max,[Ao ; Co]);
dif_canon = 1;                  % use canonical form of DE
[t_sol1, y_sol1] = ode45(@o2_homo_dif,0:t_max/200:t_max,[Ao ; Co]);
%% Plot
figure1 = figure('Name','2nd Order Rx Time course','NumberTitle','off');
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
set(leg,'Position',[0.673 0.374 0.182 0.123],'FontName','Times New Roman'...
    ,'EdgeColor',[0.95 0.95 0.95]);

% set(p2,'Visible','off');                        set(p3,'Visible','off');

% Draw mini plot
% axes2 = axes('Parent',figure1,'Position',[0.66 0.646 0.225 0.29], ...
%     'OuterPosition',[0.572 0.542 0.341 0.3774],'FontSize',8); 
% tmin = 9;                 tmax = 10;         % Range of t values for mini plot (sec) A=100,50
% tmin = 45;                 tmax = 50;         % Range of t values for mini plot (sec) A=25
% tmin = 90;                 tmax = 100;        % Range of t values for mini plot (sec) A=11,10

% idxmin = tmin/dt;         idxmax = tmax/dt;
% mp1 = plot(time(idxmin:idxmax),avg(idxmin:idxmax),'Parent',axes2,'Color','b');         hold on;
% mp1_dev1 = plot(time(idxmin:idxmax),avg(idxmin:idxmax)+sdev(idxmin:idxmax), ...
%     'Parent',axes2,'LineStyle','--','Color',[0.8 0.8 0.8]);
% mp1_dev0 = plot(time(idxmin:idxmax),avg(idxmin:idxmax)-sdev(idxmin:idxmax), ...
%     'Parent',axes2,'LineStyle','--','Color',[0.8 0.8 0.8]);     
% mp2 = plot(t_sol0,y_sol0(:,1),':g');
% mp3 = plot(t_sol1,y_sol1(:,1),':r');     
% mp4 = plot(avg_gtime',gA_all(1,:)','LineStyle',':','Color','k');          % Gillespie trajectories
                                                                                       hold off;
% set(axes2,'YMinorTick','on','XMinorTick','on','Box','off');
% xlim(axes2,[tmin tmax]);
% ylim(axes2,[0 10]);                             % for A = 100
% ylim(axes2,[2.5 7.5]);                           % for A = 50
% ylim(axes2,[0 2.5]);                             % for A = 25, 24
% ylim(axes2,[0.3 1.4]);                           % for A = 11, 10
% xlabel('t (sec)','FontSize',8);                 ylabel('N_A(t)','FontSize',8);

% Create rectangle
% annotation(figure1,'rectangle',[0.444 0.163 0.05436 0.06847],'EdgeColor','k'); % for A=100
% annotation(figure1,'rectangle',[0.444 0.163 0.05436 0.06847],'EdgeColor','k'); % for A=50

% Percent difference between ABK and DEab or DEcan at middle of mini plot
% w = 91;                         % DE index : A=100,50-> 45 , A=25 -> 81  , A=11,10->
% z = 9001;                       % ABK index: A=100,50-> 450, A=25 -> 3201, A=11,10-> 9001
% dper_ab  = 100*abs(avg(z)-y_sol0(w)) / y_sol0(w)
% dper_can = 100*abs(avg(z)-y_sol1(w)) / y_sol1(w)

%% Plot Extinction time distribution
bins = 20;

fig1 = figure;
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.

% hist(extt(:,:),bins);                         % plot histogram

[num1, c1] = hist(extt,bins);                   % histogram data 
Num1 = num1 / (reps*Ao);                        % Normalize

b = bar(c1',Num1');                               hold on;
set(b(1),'FaceColor','r','EdgeColor',[0.70 0.78 1]);        

timeD = 0:t_max/bins:t_max;

dif_canon = 0;                  % use agent-based form of DE
[t_sol2,y_sol2] = ode45(@o2_homo_dif,timeD,[Ao ; Co]);

dA = abs(diff(y_sol2(:,1))) / Ao;                     

timeDm = t_max/bins/2:t_max/bins:t_max;
plot(timeDm,dA','Color',[0 0.75 0],'LineWidth',2);

title(['N_{A,i} = ' num2str(Ao)],'FontName','Times New Roman',...
    'FontSize',12,'FontWeight','Normal');
axis([0 t_max 0 0.35]);
set(gca,'XMinorTick','on','YTickLabel',{'0','0.05','0.10','0.15','0.20','0.25','0.30','0.35'},...
    'YMinorTick','on','Box','off');
xlabel('t (sec) ');             ylabel('Fraction of 2A\rightarrowX Transitions');

% leg1 = legend();
% set(leg1,'FontName','Times New Roman','FontSize',9,...
%     'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');       

% Create textbox
annotation(fig1,'textbox',[0.02 0.94 0.05 0.07],'String',{'b)'},'FontSize',12,...
    'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');

% Coefficient of Determination R^2 for histogram fit to distribution PDF.
disp(['R^2 = ' num2str(CoefDet(Num1,dA'))]);

%% Finish
toc
clear temp* i f n t x z tmin tmax idxmin idxmax fig* p* leg* At Ct *_sol* P;
% save(['o2_A=' num2str(agents) '_extt.mat']);
