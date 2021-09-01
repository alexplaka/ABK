%   Ep <--->  E            Rates: k_f, k_r (Reverse rx activated by R)
%    |  ^                  Michaelis-Menten constants: Km_f, Km_r
%    |   \   --> X -->     k_bx and k_dx: synthesis and degradation of X
%    \   \  /    \
%    -->  R  -------->                      R: Response
%   k_b   ^      k_d
%         | k_s
%         S                            S: Signal

% Simulating mutual activation switch process for R: 
% k_b: R synthesis, 1st order process wrt Ep (but Ep is not consummed)
% k_d1: R degradation, 1st order process wrt R
% k_d2: R degradation, 2nd order process wrt R, X
% k_bx: X synthesis, 0th order, UPregulated by R
% k_dx: X degradation, 1st order wrt X
% k_s: R synthesis, 1st order process wrt S
% k_f: E synthesis, MM process
% k_r: Ep synthesis, MM process wrt R (but R is not consummed)
% Km_f, Km_r: MM constants for forward and reverse rxs

% Note: Assume number of S molecules/agents is NOT changing (eg, S is DNA molecules)

% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;      clc;     tic;
rng(0);
global agents k_b k_d1 k_d2 k_s k_bx k_dx k_f k_r Km_f Km_r S;

reps = 100;                         % Repeat simulation this many times

totalTime = 2500;                   % Simulation time (sec)
dt = 1/50;                          % Constant (fixed) time step increment (sec)
t_steps = totalTime / dt;
time = 0:dt:totalTime;

agents = 100;
k_b = 0.01;                        % 0th order R synthesis rate, UPregulated by Ep
k_d1 = 0.0;                        % 1st order degradation rate (units: 1/sec)
k_d2 = 0.0075;                     % MICROscopic 2nd order degradation rate (units: 1/sec)
k_s = 0.05;                        % 0th order R synthesis rate UPregulated by S
k_bx = 0.001;                      % 0th order synthesis of X, UPregulated by R
k_dx = 0.01;                       % 1st order degradation of X
k_f = 1;                           % basal forward rate (1st order)
k_r = 0.05;                        % basal reverse rate (1st order)
Km_f = 5;                          % MICROSCOPIC Michaelis-Menten constant for forward rx
Km_r = 10;                         % MICROSCOPIC Michaelis-Menten constant for reverse rx
S = 15;                            % Assume number of S molecules/agents is NOT changing

% Initial population sizes
Ri = 0;                 Epi = 0;                Xi = 0;

% Preallocate memory
X_all = zeros(reps,t_steps+1);         
Ep_all = zeros(reps,t_steps+1);        
R_all = zeros(reps,t_steps+1);        

%% Repeat simulation 'reps' number of times
for n=1:reps
        
    fprintf('.');
    
    % Initialize - Preallocate memory
    Tr = zeros(1,t_steps);              Tx = zeros(1,t_steps);
    Tep = zeros(1,t_steps);             Te = zeros(1,t_steps);

    P_f = zeros(1,t_steps);             P_r = zeros(1,t_steps);
    P_b = zeros(1,t_steps);             P_d2 = zeros(1,t_steps);
    P_bx = zeros(1,t_steps);

    % ******** Initial conditions - Number of R, Ep, E Agents ********
    Tr(1) = Ri;                       Tx(1) = Xi;    
    Tep(1) = Epi;                     Te(1) = agents - Epi;					
    % ****************************************************************

    tempR = zeros(1,2*agents);        tempEp = zeros(1,agents); 
    tempX = zeros(1,2*agents);
    % Put "1" where agents are "alive", then randomize the array
    for c=1:Tr(1),                  tempR(c)=1;             end
    for d=1:Tep(1),                 tempEp(d)=1;            end
    for e=1:Tx(1),                  tempX(e)=1;             end
    tempR = RandArray(tempR);                   % Randomize R array
    tempEp = RandArray(tempEp);                 % Randomize Ep array
    tempX = RandArray(tempX);                   % Randomize X array
    % Markov process, so only previous and current time steps needed --> 2 rows:    
    R = [tempR ; tempR];                        % Initialize vector storing state of R agents
    Ep = [tempEp ; tempEp];                     % Initialize vector storing state of Ep agents
    E = ~ Ep;                                   % Initialize vector storing state of E agents
    % - Ep and E are complementary (E + Ep = agents)
    X = [tempX ; tempX];                        % Initialize vector storing state of X agents
    clear c d e tempR tempEp tempX;

    % [R_ss, Ep_ss] = RateCurve;           disp(['All R_ss values = ' num2str(R_ss')]);
    % bifurc(Tr(1),Tep(1));

    % *** ABK simulation ***
    t = 1;

    P_s = k_s * S * dt;             % Probability of R synthesis process (0th order)
    P_d1 = k_d1 * dt;               % Probability of R degradation (1st order)
    P_dx = k_dx * dt;               % Probability of X degradation process (1st order wrt X)

    while t <= t_steps % && Tep(t)<agents

        P_b(t) = k_b * Tep(t) * dt;      % Probability of R synthesis process (0th order)
        P_bx(t)= k_bx * Tr(t) * dt;     % Probability of X synthesis process (0th order)
        P_f(t) = k_f / (Km_f + Tep(t)) * dt;          % wrt each Ep molecule
        P_r(t) = k_r * Tr(t) / (Km_r + Te(t)) * dt;   % wrt each E molecule
        P_d2(t) = k_d2 * Tx(t) * dt;                   % wrt each R molecule

        % Take care of 0th order processes first
        if rand < P_b(t)
            tempR = find(R(1,:)==0);    % Randomly choose R agent synthesis 
            R(2,tempR(ceil(rand * size(tempR,2)))) = 1;
        end

        if rand < P_bx(t)
            tempX = find(X(1,:)==0);    % Randomly choose X agent synthesis 
            X(2,tempX(ceil(rand * size(tempX,2)))) = 1;
        end

        if rand < P_s
            tempR = find(R(1,:)==0);    % Randomly choose R agent synthesis 
            R(2,tempR(ceil(rand * size(tempR,2)))) = 1;
        end
        % End of 0th order processes

    % The following implementation also works (treating it as a 1st order rx wrt S)     
    %     for h=1:S                               % Reaction for each S molecule/agent
    %         if rand < P_s
    %             tempR = find(R(1,:)==0);         % Randomly choose R agent which becomes "alive" 
    %             R(2,tempR(ceil(rand * size(tempR,2)))) = 1;
    %         end
    %     end

        P_d = P_d1 + P_d2(t);               % Total probability of R degradation
        tempR = find(R(1,:)==1);
        for i = 1:size(tempR,2)
            if rand < P_d                   % Degradation of R, 1st + 2nd order rx
                R(2,tempR(i)) = 0;          % R agent is degraded
            end
        end

        tempE = find(E(1,:)==1);
        for j=1:size(tempE,2)
            if rand < P_r(t)
                E(2,tempE(j)) = 0;          % Conversion of E to Ep
                Ep(2,tempE(j)) = 1;
            end
        end

        tempEp = find(Ep(1,:)==1);
        for k=1:size(tempEp,2)
            if rand < P_f(t)
                Ep(2,tempEp(k)) = 0;        % Conversion of Ep to E
                E(2,tempEp(k)) = 1;
            end
        end

        tempX = find(X(1,:)==1);
        for m = 1:size(tempX,2)
            if rand < P_dx                  % Degradation of X, 1st order rx
                X(2,tempX(m)) = 0;          % X agent is degraded
            end
        end  

        Tr(t+1) = sum(R(2,:));              Tx(t+1) = sum(X(2,:));          
        Tep(t+1) = sum(Ep(2,:));            Te(t+1) = sum(E(2,:));

        R(1,:) = R(2,:);                    X(1,:) = X(2,:);
        Ep(1,:) = Ep(2,:);                  E(1,:) = E(2,:);

        t = t + 1;
    end         % end "while" loop

    % Remove unnecessary terminal 0's from arrays
    if t < t_steps
        Tr = Tr(1:t);                               Tx = Tx(1:t);                             
        Tep = Tep(1:t);                             Te = Te(1:t);             
        P_f = P_f(1:t);                             P_r = P_r(1:t);
        P_d2 = P_d2(1:t);      P_b = P_b(1:t);      P_bx = P_bx(1:t);
    end
    
    R_all(n,:) = Tr;                Ep_all(n,:) = Tep;                X_all(n,:) = Tx;

end         % end "for n=1:reps" loop

clear R Ep E X;
%% Calculate AVERAGE + SDEV Time Course

avg_R = mean(R_all);                    sdev_R = std(R_all);            cv_R = sdev_R ./ avg_R;
avg_Ep = mean(Ep_all);                  sdev_Ep = std(Ep_all);          cv_Ep = sdev_Ep ./ avg_Ep;   
avg_X = mean(X_all);                    sdev_X = std(X_all);            cv_X = sdev_X ./ avg_X;   

%% Solve DE
if exist('t','var')==0,         finaltime = totalTime;      
else                            finaltime = (t-1) * dt;         end
[t_sol, y_sol] = ode45(@actinh2_dif,0:finaltime/500:finaltime,[Ri ; Epi ; Xi]);
%% Plot [selected] time course
trial = 3;

figure('Name','Time Course','NumberTitle','off','Position',[1 1 650 406]);           hold on;

p_Rs = plot(time(1:20:end),R_all(trial,1:20:end),'r','DisplayName','ABK N_R(t)');                                
p_Eps = plot(time(1:20:end),Ep_all(trial,1:20:end),'g','DisplayName','ABK N_{Ep}(t)');
p_Xs = plot(time(1:20:end),X_all(trial,1:20:end),'b','DisplayName','ABK N_X(t)');

p_Rd = plot(t_sol,y_sol(:,1),'--','Color',[0.75 0 1],'LineWidth',2,'DisplayName','DE N_R(t)');                                       
p_Epd = plot(t_sol,y_sol(:,2),'--','Color',[0 0.5 0],'LineWidth',2,'DisplayName','DE N_{Ep}(t)');
p_Xd = plot(t_sol,y_sol(:,3),'--','Color',[0.15 1 0.75],'LineWidth',2,'DisplayName','DE N_X(t)');                                     

axis([0 finaltime 0 agents]);                               hold off;
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('t (sec)');                 ylabel('N(t)'); 

leg0 = legend([p_Rs , p_Eps , p_Xs , p_Rd , p_Epd , p_Xd]);
set(leg0,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEastOutside');                                              

%% Plot AVERAGE Time Course

figure('Name',['Avg TC, S=' num2str(S)],...
    'NumberTitle','off','Position',[1 1 650 406]);                        hold on;

p_R = plot(time(1:20:end),avg_R(1:20:end),'r','LineWidth',2,'DisplayName','<N_R(t)>_{sim}');   
plot(time(1:20:end),avg_R(1:20:end)+sdev_R(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
plot(time(1:20:end),avg_R(1:20:end)-sdev_R(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);

p_Ep = plot(time(1:20:end),avg_Ep(1:20:end),'g','LineWidth',2,'DisplayName','<N_{Ep}(t)>_{sim}'); 
% plot(time(1:20:end),avg_Ep(1:20:end)+sdev_Ep(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
% plot(time(1:20:end),avg_Ep(1:20:end)-sdev_Ep(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);

p_X = plot(time(1:20:end),avg_X(1:20:end),'b','LineWidth',2,'DisplayName','<N_{X}(t)>_{sim}');
% plot(time(1:20:end),avg_X(1:20:end)+sdev_X(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
% plot(time(1:20:end),avg_X(1:20:end)-sdev_X(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);

p_Rd = plot(t_sol,y_sol(:,1),'--','Color',[0.75 0 1],'LineWidth',2,'DisplayName','DE N_R(t)');                                       
p_Epd = plot(t_sol,y_sol(:,2),'--','Color',[0 0.5 0],'LineWidth',2,'DisplayName','DE N_{Ep}(t)');
p_Xd = plot(t_sol,y_sol(:,3),'--','Color',[0.15 1 0.75],'LineWidth',2,'DisplayName','DE N_X(t)');                                     

axis([0 time(end) 0 agents]);                                               % axis tight;    
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');                   hold off;
xlabel('t (sec)');            ylabel('N(t)');                                        

leg1 = legend([p_R , p_Ep , p_X , p_Rd , p_Epd , p_Xd]);
set(leg1,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEastOutside');                                          

%% Plot coefficient of variation
fig2 = figure('Name','Coefficient of Variation','NumberTitle','off');                  
set(fig2,'Position',[1 1 500 450]);                                     hold on;                  
fig2.PaperUnits = 'inches';
fig2.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.

p2 = plot(time(1:20:end),cv_Ep(1:20:end),'Color',[0.2 0.75 0.5],...
    'DisplayName','ABK, Ep','LineWidth',2);               
% pt2 = plot(t_sol,1./sqrt(y_sol(:,2)),':','Color',[0.1 0.5 0.75],...
%     'DisplayName','Poisson, Xp','LineWidth',2);
pt2 = plot(time(1:20:end),1./sqrt(avg_Ep(1:20:end)),':','Color',[0.5 0.75 0.2],...
    'DisplayName','Poisson, Ep','LineWidth',2);
p3 = plot(time(1:20:end),cv_X(1:20:end),'Color',[0.35 0.5 0.75],...
    'DisplayName','ABK, X','LineWidth',2);               
pt3 = plot(time(1:100:end),1./sqrt(avg_X(1:100:end)),':','Color',[0 0.5 0.75],...
    'DisplayName','Poisson, X','LineWidth',2);
p1 = plot(time(1:20:end),cv_R(1:20:end),'k','DisplayName','ABK, R','LineWidth',2);               
pt1 = plot(time(1:20:end),1./sqrt(avg_R(1:20:end)),':','Color',[0.5 0.5 0.5],...
    'DisplayName','Poisson, R','LineWidth',2);

% set([p(2) p(3) pt(2) pt(3)],'Visible','off');       % Don't show for ...
axis([0 totalTime 0 0.6]);        
% axis tight;
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');

leg2 = legend([p1 pt1 p2 pt2 p3 pt3]);
set(leg2,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95],'Location','Best');

xlabel('t (sec)');         
ylabel('$\eta$','Interpreter','LaTeX');                                 hold off;                                    

%% Finish
fprintf(['\n' ElapsedTime(toc)]);

clear array R Ep dR dEp R_nc Ep_nc *_sym P_*;
clear h i j k n r t temp tempR tempEp tempE fig0 leg0 fig* leg* p*;
