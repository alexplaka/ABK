%    k3     k4
%    -->  X -->                        
%   /    /
%   S   /                               S: Signal
%   \  /
%    -->  R -->                         R: Response
%    k1     k2
% Simulating birth-death process for R and X, and X inhibits the production of R: 
% k1: synthesis of R, 0th order process (upregulated by S; inhibited by X)
% k2: degradation of R, 1st order process (R)
% k3: synthesis of X, 0th order process (upregulated by S)
% k4: degradation of X, 1st order process (X)

% Note: Assume number of S molecules/agents is NOT changing (eg, S is DNA molecules)
% Note: Assume number of X molecules/agents is NOT changing while inhibiting formation of R

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;     clc;     tic;
rng(0);
pb = waitbar(0,'0');

global k1 k2 k3 k4 S K nH;

k1 = 0.50;          % MICROscopic rate of synthesis of R; 0th order [units: 1/sec]
k2 = 0.10;         % MICROscopic rate of degradation of R; 2nd order [units: 1/sec]
k3 = 0.10;          % MICROscopic rate of synthesis of X; 0th order [units: 1/sec]
k4 = 0.05;          % MICROscopic rate of degradation of X; 1st order [units: 1/sec]

K = 40;             % K_50 for X inhibiting production of R (# of X to produce half-maximal effect) 
nH = 4;             % Hill coefficient for X inhibiting the production of R

agents = 200;                  % Maximum number of agents for either R or X

totalTime = 300;               % Simulation time (sec)
dt = 1/100;                    % Constant (fixed) time step increment (sec)
t_steps = totalTime / dt;
time = zeros(1,t_steps);

reps = 500;
    
% Set up step function for escalating values of S
S_level = linspace(5,50,4);
S_array = zeros(1,t_steps);
S_array(1 : t_steps/5) = S_level(1);
S_array(t_steps*1/5 + 1 : t_steps*2/5) = S_level(2);
S_array(t_steps*2/5 + 1 : t_steps*3/5) = S_level(3);
S_array(t_steps*3/5 + 1 : t_steps) = S_level(4);

% Preallocate memory for repeated simulation time courses.
R_all = zeros(reps,t_steps);                    X_all = zeros(reps,t_steps);

%% Theoretical steady-state values
R_ss = (k1 / k2) .* S_level ./ (1 + (k3 / k4 / K)^nH .* S_level.^nH);   % Steady-state value of R

X_ss(1) = k3 * S_array(t_steps/5) / k4;           % Steady-state value of X
X_ss(2) = k3 * S_array(t_steps*2/5) / k4;         % Steady-state value of X
X_ss(3) = k3 * S_array(t_steps*3/5) / k4;         % Steady-state value of X
X_ss(4) = k3 * S_array(t_steps) / k4;             % Steady-state value of X

disp(['Theoretical R_ss = ' num2str(R_ss)]);
disp(['Theoretical X_ss = ' num2str(X_ss)]);

%% ABK simulation
for n=1:reps

    progress = n/reps;
    waitbar(progress,pb,sprintf('%.0f%% %',progress*100));

    % Initialize
    Tr = zeros(1,t_steps);              Tx = zeros(1,t_steps);         
    Tr(1) = 0;                          Tx(1) = 0;              % ** Initial condition **

    tempR = zeros(1,agents);             
    % Put "1" where agents are "alive", then randomize the array
    for c=1:Tr(1),                  tempR(c)=1;             end
    tempR = RandArray(tempR);                   % Randomize array
    % Markov process, so only previous and current time steps needed --> 2 rows:    
    R = [tempR ; tempR];                        % Initialize vector storing state of R agents

    tempX = zeros(1,agents);             
    % Put "1" where agents are "alive", then randomize the array
    for d=1:Tx(1),                  tempX(d)=1;             end
    tempX = RandArray(tempX);                   % Randomize array
    % Markov process, so only previous and current time steps needed --> 2 rows:    
    X = [tempX ; tempX];                        % Initialize vector storing state of X agents
    clear c d tempR tempX;

    % Preallocate memory for Time-Dependent transition Probability
    P1 = zeros(1,t_steps-1);

    for t = 2:t_steps

%         if Tx(t-1)>=agents,         break;          end    
    %     if Tr(t-1)>=agents || Tx(t-1)>=agents,         break;          end
    %     if abs(Tr(t-1)-R_ss)<0.5 && abs(Tx(t-1)-X_ss)<0.5,      break;          end

        S = S_array(t-1);
    
%         Choose only one form of P1   !!!  0th or 1st order wrt S?
        P1(t-1) = (k1 * S) / (1 + (Tx(t-1)/K)^nH) * dt;     % P of R synthesis (0th order)
%         P1(t-1) = k1 / (1 + (Tx(t-1)/K)^nH) * dt;           % P_dif of R synthesis (1st order wrt S)
        
        P2 = 1 - exp(-k2*dt);                               % P_ber of R degradation (1st order wrt R)
    %     P2 = k2 * dt;                                       % P_dif of R degradation (1st order)
        P3 = (k3 * S) * dt;                                 % P of X synthesis (0th order)
        P4 = 1 - exp(-k4 * dt);                             % P_ber of X degradation (1st order wrt X)
    %     P4 = k4 * dt;                                       % P_dif of X degradation (1st order wrt X)

%       ** Take care of 0th order processes first **
        if rand < P1(t-1)
            temp = find(R(1,:)==0);             % Randomly choose R agent which becomes "alive" 
            R(2,temp(ceil(rand * size(temp,2)))) = 1;
        end

        if rand < P3
            temp = find(X(1,:)==0);             % Randomly choose X agent which becomes "alive" 
            X(2,temp(ceil(rand * size(temp,2)))) = 1;  
        end
%       **                                        **
        
%       Or, do reaction with k1 as 1st order wrt S (instead of 0th)
%         for z=1:S
%             if rand < P1(t-1)
%                 temp = find(R(1,:)==0);         % Randomly choose R agent which becomes "alive" 
%                 R(2,temp(ceil(rand * size(temp,2)))) = 1;
%             end
%         end
        

        tempR = find(R(1,:)==1);
        for h=1:size(tempR,2)                          % Reaction for each R molecule/agent
            if rand < P2
                R(2,tempR(h)) = 0;                     % Degradation of R
            end
        end

        tempX = find(X(1,:)==1);
        for i = 1:size(tempX,2)                        % Reaction for each X molecule/agent
            if rand < P4
                X(2,tempX(i)) = 0;                     % Degradation of X                         
            end
        end

        Tr(t) = sum(R(2,:));                Tx(t) = sum(X(2,:));
        R(1,:) = R(2,:);                    X(1,:) = X(2,:);
        time(t) = time(t-1) + dt;

    end         % end "for t=2:t_steps" loop
    
    R_all(n,:) = Tr;                    X_all(n,:) = Tx;

end             % end "for n=1:reps" loop

avgR = mean(R_all);                  sdevR = std(R_all);
avgX = mean(X_all);                  sdevX = std(X_all);

% Remove unnecessary terminal 0's from arrays
% if t < t_steps
%     Tr = Tr(1:t-1);             Tx = Tx(1:t-1);                 P2 = P2(1:t-1);  
%     time = time(1:t-1);         S_array = S_array(1:t-1);                     
% end

disp(['ABK sim terminal <R> value  = ' num2str(avgR(end))]);  
% disp(['ABK sim terminal <X> value  = ' num2str(avgX(end))]);        

tmax = time(end);
%% Solve DE
% [t_sol, y_sol] = ode45(@adapt_dif,0:tmax/500:tmax,[Tr(1) Tx(1)]);     % General - constant S

S = S_level(1);                                                         % First Level of S
[t_sol1, y_sol1] = ode45(@adapt_dif,0:totalTime/5,[Tr(1) Tx(1)]);

S = S_level(2);                                                         % Second Level of S
[t_sol2, y_sol2] = ode45(@adapt_dif,totalTime/5+dt:totalTime*2/5,[y_sol1(end,1) y_sol1(end,2)]);

S = S_level(3);                                                         % Third level of S
[t_sol3, y_sol3] = ode45(@adapt_dif,totalTime*2/5+dt:totalTime*3/5,[y_sol2(end,1) y_sol2(end,2)]);

S = S_level(4);                                                         % Fourth Level of S
[t_sol4, y_sol4] = ode45(@adapt_dif,totalTime*3/5+dt:totalTime,[y_sol3(end,1) y_sol3(end,2)]);

t_sol = [t_sol1; t_sol2; t_sol3; t_sol4];
y_sol = [y_sol1; y_sol2; y_sol3; y_sol4];

%% Plot Time Course
fig1 = figure('Name','Time Course - Adaptation','NumberTitle','off');
set(fig1,'Position',[1 1 500 406]);

p1 = plot(time,avgR,'b','MarkerSize',1,'DisplayName','<N_R(t)>_{sim}');                 hold on;
% p1_dev1 = plot(time,avgR+sdevR,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p1_dev0 = plot(time,avgR-sdevR,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p1s = plot(time,R_all(end-4,:),'m','MarkerSize',3,'DisplayName','Sample Trajectory');

p2 = plot(time,avgX,':','Color',[0 0.8 0],'MarkerSize',3,'DisplayName','<N_X(t)>_{sim}');                 
% p2_dev1 = plot(time,avgX+sdevX,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p2_dev0 = plot(time,avgX-sdevX,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p2s = plot(time,X_all(end,:),'m','MarkerSize',3,'DisplayName','Sample Trajectory');

p1d = plot(t_sol,y_sol(:,1),'--c','LineWidth',1,'DisplayName','DE N_R(t)');
% p2d = plot(t_sol,y_sol(:,2),'--g','LineWidth',1,'DisplayName','DE N_X');

pS = plot(time,S_array,'r','LineWidth',2,'DisplayName','N_S');

for i=1:size(R_ss,2)-1
    time_part = totalTime*(i-1)/5:totalTime*i/5;
    p_ss(i) = plot(time_part,R_ss(i)*ones(1,size(time_part,2)),'--',...
        'Color',[1 0.75 0.2],'LineWidth',2,'DisplayName','N_R^*');
end
time_part = totalTime*(i)/5:totalTime;   % Plot R_ss for last value of S - Extended time interval
p_ss(i+1) = plot(time_part,R_ss(i+1)*ones(1,size(time_part,2)),'--',...
    'Color',[1 0.75 0.2],'LineWidth',2,'DisplayName','N_R^*');

axis([0 tmax 0 100]);         
xlabel('t (sec)');        ylabel('N(t)');                                           hold off;
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');

% leg1 = legend('R stoc','X stoc','R deter','X deter','S','R_{ss}');
leg1 = legend([p1 p2 p1d pS p_ss(1)]);
set(leg1,'Location','NorthWest');
set(leg1,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95]);

%% Finish
close(pb);
clear pb h i n t r w p* temp tempX tempR X R;
toc
