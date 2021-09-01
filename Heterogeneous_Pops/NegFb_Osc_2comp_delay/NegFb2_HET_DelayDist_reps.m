%      R <---->  Rp             Rates: k_f, k_r (Forward rx activated by Yp)
%           ^     |             Michaelis-Menten constants: Km_f, Km_r
%          /      |
%         /       |
%        /        |
%        |        \
%   -->  X  ----------->                      Rp: Response
%   k_b  ^    k_d1/k_d2
%        | 
%        | k_s
%        S                             S: Signal

% Simulating negative feedback process (2-component system: X, R/Rp) 
% k_b: X synthesis, 0th order process
% k_s: X synthesis, 1st order process wrt S (but S is not consummed)
% k_s: ALTERNATIVE IMPLEMENTATION: X synthesis, 0th order process, UPregulated by S
% k_d1: X degradation, 1st order process wrt X
% k_d2: X degradation, 2nd order process wrt X, Rp (but Rp is not consummed)
% k_f: Rp synthesis, MM process wrt Yp (but Yp is not consummed)
% k_r: R synthesis, MM process 
% Km_f, Km_r: MM constants for forward and reverse rxs R <==> Rp

% ** We simulate a system with delay in the negative feedback rx only. **

% Note: Assume number of S molecules/agents is NOT changing (eg, S is DNA molecules)

% Here, we consider a normal distribution of delay values.
% All X agents have a delay value associated with them.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;       
clc;     tic;                                          
rng(1);
%% Declare Parameters
global agents k_s k_b k_d1 k_d2 k_f k_r Km_f Km_r S;

reps = 100;                             % Repeat simulation this many times

totalTime = 2000;                       % Simulation time (sec)
dt = 0.01;                              % Constant (fixed) time step increment (sec)
t_steps = totalTime / dt;
time = 0:dt:totalTime;
to = 100;

lag_mean = 50;
lag_std = 10;

agents = 100;
k_b = 0;                            % 0th order X synthesis rate
k_d1 = 0;                           % 1st order degradation rate (units: 1/sec)
k_d2 = 0.005;                       % 2nd order degradation rate wrt X, Rp (units: 1/sec)
k_s = 0.10;                         % 0th order R synthesis rate UPregulated by S
k_f = 0.01;                         % basal forward rate (1st order)
k_r = 1;                            % basal reverse rate (1st order)
Km_f = 10;                          % MICROSCOPIC Michaelis-Menten constant for forward rx
Km_r = 10;                          % MICROSCOPIC Michaelis-Menten constant for reverse rx

S = 30;                             % Assume number of S molecules/agents is NOT changing

% Initial population sizes
Xi = 0;
Rpi= 0;                         Ri = agents - Rpi;

% Preallocate memory
X_all = zeros(reps,t_steps+1);         
Rp_all = zeros(reps,t_steps+1);        

%% Repeat simulation
for n=1:reps
        
    fprintf('.');
    
    % Initialize, set Initial conditions.
    Tx = zeros(1,t_steps);
    Tr = zeros(1,t_steps);              Trp = zeros(1,t_steps);
    P_f = zeros(1,t_steps);             P_r = zeros(1,t_steps);
    P_d2 = zeros(1,t_steps);

    % ******** Initial conditions - Number of R, Ep, E Agents ********
    Tx(1) = Xi;    
    Trp(1) = Rpi;                       Tr(1) = agents - Trp(1);					
    % ****************************************************************

    tempX = zeros(1,3*agents);          tempR = zeros(1,agents);             
    % Put "1" where agents are "alive", then randomize the array
    for c=1:Tx(1),                      tempX(c)=1;             end
    for e=1:Tr(1),                      tempR(e)=1;             end
    tempX = RandArray(tempX);                    % Randomize X array
    tempR = RandArray(tempR);                    % Randomize R array
%   Markov process, so only previous and current time steps needed --> 2 rows:    
    Xv = [tempX ; tempX];                        % Initialize Vector storing state of X agents
    Rv = [tempR ; tempR];                        % Initialize Vector storing state of R agents
    Rpv = ~ Rv;                                  % Initialize Vector storing state of Rp agents
    % - R and Rp are complementary (R + Rp = agents).
    clear c d e tempX tempR;
        
% ABK simulation
    t = 1;

    P_b = k_b * dt;                  % Probability of X synthesis process (0th order)
    P_s = k_s * S * dt;              % Probability of X synthesis process (0th order), UPreg'd by S
    P_d1 = 1 - exp(-k_d1 * dt);      % Probability of X degradation process (1st order) wrt X

    while t <= t_steps % && Tep(t)<agents

        % ****** Introducing Delay in the Probability Expression of Feedback Rx ****** 
        if t <= to/dt
            P_d2(t) = k_d2 * Trp(t) * dt;                % P wrt each X molecule
        else    % introduce normally-distributed delay of lag_mean seconds, with sdev = lag_std
            P_d2(t) = k_d2 * Trp(t - round(lag_mean + lag_std*randn) /dt) * dt;           
        end                        % round to nearest integer # of seconds
        % ****************************************************************************
        
        P_f(t)  = k_f * Tx(t) / (Km_f + Tr(t)) * dt;          % wrt each R molecule
        P_r(t)  = k_r / (Km_r + Trp(t)) * dt;                 % wrt each Rp molecule

        % Take care of 0th order processes first
        tempXd = find(Xv(1,:)==0);                
        if rand < P_b            
            chooseX1 = ceil(rand * size(tempXd,2));         % Randomly choose X agent synthesis 
            Xv(2,tempXd(chooseX1)) = 1;  
            
            % Now, ensure newly created agent is not re-created in the same time step
            % by the other 0th order process.
            tempXd(chooseX1) = [];                          
        end

        if rand < P_s
            chooseX2 = ceil(rand * size(tempXd,2));         % Randomly choose X agent synthesis 
            Xv(2,tempXd(chooseX2)) = 1;
        end
        % End of 0th order processes

    % The following implementation also works (treating it as a 1st order rx wrt S)     
    %     for h=1:S                               % Reaction for each S molecule/agent
    %         if rand < P_s
    %             tempX = find(Xv(1,:)==0);       % Randomly choose X agent which becomes "alive" 
    %             Xv(2,tempR(ceil(rand * size(tempX,2)))) = 1;
    %         end
    %     end

        P_d = P_d1 + P_d2(t);               % Total probability of X degradation
        
        tempX = find(Xv(1,:)==1);        
        for i = 1:size(tempX,2)
            if rand < P_d                   % Degradation of X, 2 contributing processes
                Xv(2,tempX(i)) = 0;         % X agent is degraded
            end
        end

        tempR = find(Rv(1,:)==1);
        for j=1:size(tempR,2)
            if rand < P_f(t)
                Rv(2,tempR(j)) = 0;         % Conversion of R to Rp
                Rpv(2,tempR(j)) = 1;
            end
        end

        tempRp = find(Rpv(1,:)==1);
        for k=1:size(tempRp,2)
            if rand < P_r(t)
                Rpv(2,tempRp(k)) = 0;       % Conversion of Rp to R
                Rv(2,tempRp(k)) = 1;
            end
        end

        Tx(t+1) = sum(Xv(2,:));
        Tr(t+1) = sum(Rv(2,:));             Trp(t+1) = sum(Rpv(2,:));      

        Xv(1,:) = Xv(2,:);
        Rv(1,:) = Rv(2,:);                  Rpv(1,:) = Rpv(2,:);
        t = t + 1;
    end

    X_all(n,:) = Tx;                      Rp_all(n,:) = Trp;
        
end             % *** end 'for n=1:reps' loop ***

clear Tx Tr Trp Xv Rv Rpv temp* i j k n t chooseX*;
clear P_b P_s P_d P_d1 P_d2 P_f P_r;

%% Calculate AVERAGE + SDEV Time Course  ** (Start from here when loading data file) **

avg_X = mean(X_all);                    sdev_X = std(X_all);            cv_X = sdev_X ./ avg_X;   
avg_Rp = mean(Rp_all);                  sdev_Rp = std(Rp_all);          cv_Rp = sdev_Rp ./ avg_Rp;

N_max = max([max(avg_X+sdev_X) max(avg_Rp+sdev_Rp)]);

%% Solve Delay DE 
global DEhistory;
% Solve ODE for the initial time period 0:to
t_sol0 = [0:to]';
DEhistory = ode23(@negFb_dif,t_sol0,[Xi; Rpi]);   % DEhistory is a structure

y_sol0 = deval(DEhistory,t_sol0)';

cd DDE_files;

DDEsol = dde23(@negFb_delaydif,lag_mean,@negFb2_dde_hist,[to , totalTime]);

% plot(DDEsol.x,DDEsol.y);

t_sol1 = [to:totalTime]';
y_sol1 = deval(DDEsol,t_sol1)';                           cd ..;

t_sol = [t_sol0 ; t_sol1];
y_sol = [y_sol0 ; y_sol1];

%% Plot [selected] time course
% Using Report Nomenclature: X = R , Rp = Xp  (script name = report name) in legend.

trial = 1;                         % Choose simulation run/trial (trial <= reps)

fig0 = figure('Name','Time Course','NumberTitle','off');                    hold on;
set(fig0,'Position',[1 1 500 450]);
set(fig0,'PaperPosition',[0 0 6 5],'PaperUnits','inches');

plot(time(1:20:end),X_all(trial,1:20:end),'b'); 
plot(time(1:20:end),Rp_all(trial,1:20:end),'r');
plot(t_sol,y_sol(:,1),'--','Color',[0.15 1 0.75],'LineWidth',2); 
plot(t_sol,y_sol(:,2),'--','Color',[0.75 0 1],'LineWidth',2);  

% axis([0 time(end) 0 N_max]);  
axis([0 time(end) 0 165]);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('$t \; \textrm{(sec)}$','Interpreter','LaTeX');                 
ylabel('$N(t)$','Interpreter','LaTeX');                                       hold off;

title(['$\textbf{Normally-distributed lag: } \bar{\delta} =' num2str(lag_mean) ...
    '\, \textrm{sec} \, \textrm{,} \; \sigma_{\delta} =' num2str(lag_std) ...
    '\, \textrm{sec} $'],'Interpreter','LaTeX');

leg0 = legend('$\textrm{ABK } N_R(t)$','$\textrm{ABK } N_{Xp}(t)$',...
    '$\textrm{DDE } N_R(t)$','$\textrm{DDE } N_{Xp}(t)$','Location','NorthEast');
% leg0 = legend('ABK Sample Run N_R(t)','ABK Sample Run N_{Xp}(t)','Location','Best'); % Without DEs
set(leg0,'Interpreter','LaTeX','FontSize',10,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');

% Create textbox
annotation(fig0,'textbox',...
    [0.00771 0.93471 0.05010 0.06429],'String','a)',...
    'FontSize',12,'FontName','Times New Roman',...
    'FitBoxToText','off','EdgeColor','none');

%% Plot AVERAGE Time Course
% Using Report Nomenclature: X = R , Rp = Xp  (script name = report name) in legend.

fig1 = figure('Name','Avg TC','NumberTitle','off');                              hold on;
set(fig1,'Position',[1 1 500 450]);
set(fig1,'PaperPosition',[0 0 6 5],'PaperUnits','inches');

p_avgX = plot(time(1:20:end),avg_X(1:20:end),'b','LineWidth',2,...
    'DisplayName','$\textrm{ABK } < \! N_R(t) \! >$');
plot(time(1:20:end),avg_X(1:20:end)+sdev_X(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
plot(time(1:20:end),avg_X(1:20:end)-sdev_X(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);

p_avgRp = plot(time(1:20:end),avg_Rp(1:20:end),'r','LineWidth',2,...
    'DisplayName','$\textrm{ABK } < \! N_{Xp}(t) \! >$');
% plot(time(1:20:end),avg_Rp(1:20:end)+sdev_Rp(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
% plot(time(1:20:end),avg_Rp(1:20:end)-sdev_Rp(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);

p_ddeX = plot(t_sol,y_sol(:,1),'--','Color',[0.15 1 0.75],...
    'LineWidth',2,'DisplayName','$\textrm{DDE } N_R(t)$'); 
p_ddeRp = plot(t_sol,y_sol(:,2),'--','Color',[0.75 0 1],...
    'LineWidth',2,'DisplayName','$\textrm{DDE } N_{Xp}(t)$');  

% axis([0 time(end) 0 1.05*N_max]);                                           % axis tight;   
axis([0 time(end) 0 165]); 
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('$t \; \textrm{(sec)}$','Interpreter','LaTeX');                 
ylabel('$N(t)$','Interpreter','LaTeX');                                       hold off;

title(['$\textbf{Normally-distributed lag: } \bar{\delta} =' num2str(lag_mean) ...
    '\, \textrm{sec} \, \textrm{,} \; \sigma_{\delta} =' num2str(lag_std) ...
    '\, \textrm{sec} $'],'Interpreter','LaTeX');

leg1 = legend([p_avgX , p_avgRp , p_ddeX , p_ddeRp]);
% leg1 = legend([p_avgX , p_avgRp]);            % without DEs
set(leg1,'Interpreter','LaTeX','FontSize',10,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');                       hold off;

% Create textbox
annotation(fig1,'textbox',...
    [0.00771 0.93471 0.05010 0.06429],'String','b)',...
    'FontSize',12,'FontName','Times New Roman',...
    'FitBoxToText','off','EdgeColor','none');

%% Plot coefficient of variation
fig2 = figure('Name','Coefficient of Variation','NumberTitle','off');                  
set(fig2,'Position',[1 1 500 450]);                                     hold on;                  
fig2.PaperUnits = 'inches';
fig2.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.

p1 = plot(time,cv_X,'k','DisplayName','ABK R','LineWidth',2);               
pt1 = plot(time,1./sqrt(avg_X),':','Color',[0.2 0.75 0.5],...
    'DisplayName','Poisson, R','LineWidth',2);
p2 = plot(time,cv_Rp,'Color',[0.5 0.5 0.5],'DisplayName','ABK Xp','LineWidth',2);               
pt2 = plot(time,1./sqrt(avg_Rp),':','Color',[0.1 0.5 0.75],...
    'DisplayName','Poisson, Xp','LineWidth',2);

% set([p(2) p(3) pt(2) pt(3)],'Visible','off');       % Don't show for ...
axis([0 totalTime 0 1]);        
% axis tight;
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');

leg2 = legend([p1 pt1 p2 pt2]);
set(leg2,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95],'Location','Best');

xlabel('t (sec)');         
ylabel('$\eta$','Interpreter','LaTeX');                                 hold off;                                    

% savefig('CV.fig');

%% Finish
% save('matlab.mat','k*','K*','S','agents','lag','Xi','Rpi','reps','totalTime','dt','time','*_all');

fprintf(['\n' ElapsedTime(toc)]);

clear array X Rp dX dRp X_nc Rp_nc X_nc_sym Rp_nc_sym Yp leg* fig* T* P* p* Yv Ypv Rv Rpv Xv;
clear m h i j k n q r t J Jac lambdas temp tempR tempRp tempX tempY tempYp *_sym ans;  
