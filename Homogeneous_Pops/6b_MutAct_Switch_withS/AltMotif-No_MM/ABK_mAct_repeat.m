%   Ep  <--->  E            Rates: k_f, k_r (Reverse rx activated by R)
%    |   ^                  
%    |    \
%    \    |
%    -->  R   -->                      R: Response
%   k_b   ^   k_d
%         | k_s
%         S                            S: Signal

% Simulating mutual activation switch process for R: 
% k_b: R synthesis, 1st order process wrt Ep (but Ep is not consummed)
% k_b: ALTERNATIVE IMPLEMENTATION: R synthesis, 0th order process, UPregulated by Ep
% k_d: R degradation, 1st order process wrt R
% k_s: R synthesis, 1st order process wrt S
% k_s: ALTERNATIVE IMPLEMENTATION: R synthesis, 0th order process, UPregulated by S
% k_f: E synthesis
% k_r: Ep synthesis (but R is not consummed)

% Note: Assume number of S molecules/agents is NOT changing (eg, S is DNA molecules)

% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear;      clc;     tic;                           
rng(1);
reps = 30;                 % Number of times to repeat experiment
global agents k_b k_d k_s k_f k_r S;

totalTime = 10;              % Simulation time (sec)
dt = 0.01;                   % Constant (fixed) time step increment (sec)
t_steps = totalTime / dt;
time = zeros(1,t_steps);

agents = 100;
k_b = 0.01;                        % 0th order R synthesis rate, UPregulated by Ep
k_d = 0.1;                        % 1st order degradation rate (units: 1/sec)
k_s = 0.05;                        % 0th order R synthesis rate UPregulated by S
k_f = 0.2;                         % basal forward rate (1st order)
k_r = 0.05;                        % basal reverse rate (1st order)
S = 10;                            % Assume number of S molecules/agents is NOT changing

Ri = 60;            Epi = 30;      % Initial conditions

% Clumped parameters for determining R_ss
a = 1;
b = - k_s/k_d * S + k_f/k_r - k_b/k_d * agents;
c = - k_f*k_s / (k_r*k_d);
R_ss_exact(1) = (-b + sqrt(b^2 - 4*a*c)) / (2*a);
R_ss_exact(2) = (-b - sqrt(b^2 - 4*a*c)) / (2*a);
R_ss = R_ss_exact(find(R_ss_exact>0));
disp(['Theory: R_ss = ' num2str(R_ss)]);
Ep_ss = (k_d * R_ss - k_s * S) / k_b;

for n = 1:reps
    
    fprintf(1,'.');
    
    Tr = zeros(1,t_steps);
    Tep = zeros(1,t_steps);             Te = zeros(1,t_steps);
    
    P_f = zeros(1,t_steps);             P_r = zeros(1,t_steps);
    P_b = zeros(1,t_steps);

    % ******** Initial conditions - Number of R, Ep, E Agents ********
    Tr(1) = Ri;        Tep(1) = Epi;      Te(1) = agents - Tep(1);					
    % ****************************************************************

    tempR = zeros(1,agents);        tempEp = zeros(1,agents);             
    % Put "1" where agents are "alive", then randomize the array
    for c=1:Tr(1),                  tempR(c)=1;             end
    for d=1:Tep(1),                 tempEp(d)=1;            end
    tempR = RandArray(tempR);                   % Randomize R array
    tempEp = RandArray(tempEp);                 % Randomize Ep array
    % Markov process, so only previous and current time steps needed --> 2 rows:    
    R = [tempR ; tempR];                        % Initialize vector storing state of R agents
    Ep = [tempEp ; tempEp];                     % Initialize vector storing state of Ep agents
    E = ~ Ep;                                   % Initialize vector storing state of E agents
    % - Ep and E are complementary (E + Ep = agents)
    clear c d tempR tempEp;

    t = 1;

    P_d = k_d * dt;                 % Probability of R degradation process (1st order wrt R)
    P_s = k_s *S* dt;               % Probability of R synthesis process (0th order)

    while t <= t_steps % && Tep(t)<agents
        P_b(t) = k_b * Tep(t) * dt;      % Probability of R synthesis process (0th order)
        P_f(t) = k_f * dt;              % wrt each Ep molecule
        P_r(t) = k_r * Tr(t) * dt;   % wrt each E molecule

        % Take care of 0th order processes first
        if rand < P_b(t)
            tempR = find(R(1,:)==0);    % Randomly choose R agent synthesis 
            R(2,tempR(ceil(rand * size(tempR,2)))) = 1;
        end

        if rand < P_s
            tempR = find(R(1,:)==0);    % Randomly choose R agent synthesis 
            R(2,tempR(ceil(rand * size(tempR,2)))) = 1;
        end
        % End of 0th order processes

        tempR = find(R(1,:)==1);
        for i = 1:size(tempR,2)
            if rand < P_d                   % Degradation of R, 1st order rx
                R(2,tempR(i)) = 0;          % R agent is degraded
            end
        end

        tempE = find(E(1,:)==1);
        for j=1:size(tempE,2)
            if rand < P_r(t)
                E(2,tempE(j)) = 0;          % Conversion of E to Ep
                tempEp = find(Ep(1,:)==0);  % Randomly choose Ep agent synthesis
                Ep(2,tempEp(ceil(rand * size(tempEp,2)))) = 1;
            end
        end

        tempEp = find(Ep(1,:)==1);
        for k=1:size(tempEp,2)
            if rand < P_f(t)
                Ep(2,tempEp(k)) = 0;        % Conversion of Ep to E
                tempE = find(E(1,:)==0);    % Randomly choose E agent synthesis
                E(2,tempE(ceil(rand * size(tempE,2)))) = 1;
            end
        end

        Tr(t+1) = sum(R(2,:));          
        Tep(t+1) = sum(Ep(2,:));            Te(t+1) = sum(E(2,:));
        
        R(1,:) = R(2,:);
        Ep(1,:) = Ep(2,:);                  E(1,:) = E(2,:);
        
        t = t + 1;
    end
    % Remove unnecessary terminal 0's from arrays
    if t < t_steps
        Tr = Tr(1:t);                                         
        Tep = Tep(1:t);                             Te = Te(1:t);             
        P_f = P_f(1:t);                             P_r = P_r(1:t);                           
    end

    Rv(n,:) = Tr;                           Epv(n,:) = Tep;

    finaltime = (t-1) * dt;

end                 % end 'for n'

[t_sol, y_sol] = ode45(@mAct_dif,0:finaltime/500:finaltime,[Ri; Epi; agents-Epi]);
%% Plot (selected) time trajectories
time = 0:dt:finaltime;
scatter(t_sol,y_sol(:,1),3,'xr');                                       hold on;
% scatter(t_sol,y_sol(:,2),3,'.r');
% scatter(t_sol,y_sol(:,3),3,'.g'); 
scatter(time,Rv(7,:),3,'.b');                                                
% scatter(time,Tep,3,'or');
% scatter(time,Te,3,'og');
axis([0 finaltime 0 agents]);         xlabel('time');            ylabel('#');     hold off; 
% legend('R stoc','Ep stoc','E stoc','R deter','Ep deter','E deter','Location','Best');

%% Calculate AVERAGE + SDEV Time Course
avg_R = mean(Rv);                       sdev_R = std(Rv);
avg_Ep = mean(Epv);                     sdev_Ep = std(Epv);
%% Plot AVERAGE Time Course - SINGLE PLOT

fig2 = figure('Name','Avg TC','NumberTitle','off','Position',[1 1 500 500]);  
title(['N_{R,i} = ' num2str(Ri) ' , N_{Ep,i} = ' num2str(Epi)],...
    'FontSize',12,'FontName','Times New Roman');                            hold on;

p_R = plot(time(1:20:end),avg_R(1:20:end),'r','LineWidth',3,'DisplayName','<N_{R}(t)>_{sim}');   
plot(time(1:20:end),avg_R(1:20:end)+sdev_R(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
plot(time(1:20:end),avg_R(1:20:end)-sdev_R(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);

p_Ep = plot(time(1:20:end),avg_Ep(1:20:end),'g','LineWidth',3,'DisplayName','<N_{F}(t)>_{sim}');
plot(time(1:20:end),avg_Ep(1:20:end)+sdev_Ep(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
plot(time(1:20:end),avg_Ep(1:20:end)-sdev_Ep(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);

p_Rd = plot(t_sol,y_sol(:,1),'--','Color',[0.75 0 1],'LineWidth',2,'DisplayName','DE N_R(t)');
p_Epd = plot(t_sol,y_sol(:,2),'--','Color',[0 0.5 0],'LineWidth',2,'DisplayName','DE N_{Ep}(t)');                                     

% axis tight;                                         
axis([0 time(end) 0 agents]);                                                   
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');                   hold off;
xlabel('t (sec)');            ylabel('N(t)');                                        

leg2 = legend([p_R , p_Ep , p_Rd , p_Epd]);
set(leg2,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');                                          

%% Plot State Space - Modify accordingly
figure('Name','State Space','NumberTitle','off');               hold on;
plot(y_sol(:,1),y_sol(:,2),'+k');                 % plot deterministic trajectory

% Modify the following for plotting specific trajectories
scatter(Rv(29,:),Epv(29,:),3,'.r');        % plot sample stochastic trajectory
scatter(avg_R,avg_Ep,3,'.b');              % plot stochastic trajectory

plot(R_ss,Ep_ss,'oc','MarkerSize',6);
xlabel('R');                    ylabel('E_p');                  hold off; 
legend('Deterministic','ABM stochastic');

%% Finish
clear h i j k n t Tr Tep Te R Ep E temp tempR tempEp tempE;
toc