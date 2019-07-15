%    k_da           k_db
%    <--  A _  _ B  -->                 % For Report: B = R
%         ^  \/  ^
%     k_a |  /\  | k_b
%         | / ~\ |
%         |~  | ~|
%             |
%             S
% A and B are synthesized (0th order constants k_a, k_b respectively)
% A and B are degraded (k_da, k_db)
% A and B inhibit each other's synthesis (MUTUAL inhibition, alpha < 1):
% B influences the rate of the synthesis of A
% A influences the rate of the synthesis of B
% S inhibits the inhibition of B by A (assume complete repression)
% Assume S does not change.

% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

%% Declare variables and functions

clear;            tic;                  clc; 
pb = waitbar(0,'0');                      % Slows things down quite a bit

reps = 100;

global k_da k_db agents;
global k_a alpha_a K_b n_b;
global k_b alpha_b K_a n_a;
global K_s n_s S;

maxTime = 5000;                        % Maximum Simulation time (sec)
dt = 1/50;                             % Constant (fixed) time step increment (sec)
t_steps = maxTime / dt;

k_a = 1;                            % MICROSCOPIC basal rate of A synthesis (0th order)
k_b = 1;                            % MICROSCOPIC basal rate of B synthesis (0th order)
k_da = 0.01;                        % degradation rate constant (1st order) for A
k_db = 0.01;                        % degradation rate constant (1st order) for B
% Populations producing half-maximal regulatory effect
K_b = 50;                       K_a =  50;
K_s = 50;

S_array = 0:20;                      % Number of S molecules to do simulation for

% *** Feedback parameters ***
% Degree of activation: >1 activator, <1 repressor, =1 no regulation
alpha_a = 0;                        % Effect of A on synthesis of B; = 0 means Complete Repression
alpha_b = 0;                        % Effect of B on synthesis of A; = 0 means Complete Repression
% Hill Coefficient: measure of cooperativity
n_b = 3;                      n_a = 3;
n_s = 1;

% Specify initial population sizes here
Ao = 10;                       Bo = 10;

ABK_B = zeros(reps,size(S_array,2));              % Preallocate memory

%% Simulation
for a=1:reps

    progress = a / reps;
    waitbar(progress,pb,sprintf('%.0f%% %',progress*100));
    fprintf(1,'\n');

    for s=1:size(S_array,2)

        fprintf(1,'.');

        S = S_array(s);

        % Initialize - Preallocate memory for variables; ** Initial conditions **
        t = 1;                              % Time counter variable
        P_sa = zeros(1,t_steps);            % For storing probability of A synthesis at each time step
        P_sb = zeros(1,t_steps);            % For storing probability of B synthesis at each time step

        Ta = zeros(1,t_steps);              % Sum of A agents in each time step	
        Tb = zeros(1,t_steps);              % Sum of B agents in each time step

        % ******** Initial conditions - Number of A, B Agents ********
        agents = 200;                       % disp(['Agents = ' num2str(agents)]);
        Ta(1) = Ao;                         Tb(1) = Bo;			
        % ************************************************************

        tempA = zeros(1,agents);            tempB = zeros(1,agents);              
        % Put "1" where agenst are "alive", then randomize the array
        for c=1:Ta(1),                      tempA(c)=1;             end
        for d=1:Tb(1),                      tempB(d)=1;             end
        tempA = RandArray(tempA);           % Randomize A array
        tempB = RandArray(tempB);           % Randomize B array
        Av = [tempA ; tempA];               % Initialize vector storing state of A agents  
        Bv = [tempB ; tempB];               % Initialize vector storing state of B agents  
        % Notes on Av, Bv:
        % - Markov process, so only previous and current time steps needed --> 2 rows:          
        clear c d tempA tempB;

        % ABK Simulation
        P_da = 1 - exp(-k_da * dt);                % P_ber of degradation of A, 1st order wrt A
        P_db = 1 - exp(-k_db * dt);                % P_ber of degradation of B, 1st order wrt B

        while t <= t_steps % && Ta(t) > agents/100
        %     fprintf(1,'.');

            k_sa = k_a * (K_b^n_b + alpha_b * Tb(t)^n_b) / (K_b^n_b + Tb(t)^n_b);
            k_sb = k_b * (K_a^n_a + alpha_a * (Ta(t) * 1/(1+(S/K_s)^n_s))^n_a) / (K_a^n_a + (Ta(t) * 1/(1+(S/K_s)^n_s))^n_a);

            P_sa(t) = k_sa * dt;                % 0th order synthesis of A (regulated by B)
            P_sb(t) = k_sb * dt;                % 0th order synthesis of B (regulated by A)

            % Take care of 0th order processes first
            if rand < P_sa(t)
                tempA = find(Av(1,:)==0);        % Randomly choose A agent synthesis 
                Av(2,tempA(ceil(rand * size(tempA,2)))) = 1;
            end

            if rand < P_sb(t)
                tempB = find(Bv(1,:)==0);        % Randomly choose B agent synthesis 
                Bv(2,tempB(ceil(rand * size(tempB,2)))) = 1;
            end
            % end of 0th order processes

            tempA = find(Av(1,:)==1);
            for i = 1:size(tempA,2)
                if rand < P_da                   % Degradation of A, 1st order rx
                    Av(2,tempA(i)) = 0;          % A agent is degraded
                end
            end

            tempB = find(Bv(1,:)==1);
            for j = 1:size(tempB,2)
                if rand < P_db                   % Degradation of B, 1st order rx
                    Bv(2,tempB(j)) = 0;          % A agent is degraded
                end
            end

            Ta(t+1) = sum(Av(2,:));                  Tb(t+1) = sum(Bv(2,:));
            Av(1,:) = Av(2,:);                       Bv(1,:) = Bv(2,:); 
            t = t + 1;   
        end

        % Remove unnecessary terminal 0's from arrays
        % if t < t_steps
        %     Ta = Ta(1:t);                            Tb = Tb(1:t);          
        %     P_sa = P_sa(1:t);                        P_sb = P_sb(1:t);
        %     P_da = P_da(1:t);                        P_db = P_db(1:t);
        % end
        % finaltime = (t-1) * dt;

        %   Average Rp population over the last 1000 sec of the simulation
        %   Steady-state has been reached for this time period
            ABK_B(a,s) = mean(Tb(end-1000/dt:end));

    end             % end "for s=1:size(S_array,2)" loop
end                 % end "for a=1:reps" loop

avg_Bss = mean(ABK_B);
sdev_Bss = std(ABK_B);

%% Plot signal-response curve
fig1 = figure('Name','ABK SR Curve','NumberTitle','off');
set(fig1,'Position',[1 1 500 406]);                                 hold on;

% Plot Theoretical SR curve
SR = load('./Theory/SR_curve.mat','S_array','B_ss');    % Make sure file has been generated.
for h=1:size(SR.B_ss,2)
    if mod(h,2) == 0
        p(h) = plot(SR.S_array,SR.B_ss(:,h),'--',...
            'LineWidth',2,'Color',[1 0.69 0.39],'DisplayName','Unstable');
    else
    p(h) = plot(SR.S_array,SR.B_ss(:,h),'-b','LineWidth',2,'DisplayName','Stable');
    end
end

% Plot ABK results - Choose one of the following three plots
% p_ABK = plot(S_array,avg_Bss,'or','DisplayName','ABK');
p_ABK = errorbar(S_array,avg_Bss,sdev_Bss,'or',...
    'DisplayName','<N_R^*>_{sim}');                             % Using Matlab built-in function
% p_ABK = shadedErrorBar(S_array,avg_Bss,sdev_Bss,'r',1);       % Using File Exchange script

xlabel('N_S');                               
ylabel('N_R^*');                            % For Report   

% axis([0 S_array(end) 0 max(max(ABK_B))]);   
axis([S_array(1) S_array(end) 0 120]);                            
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');           hold off;

% leg1 = legend([p_ABK],'Location','SouthEast');
% set(leg1,'FontName','Times New Roman','FontSize',9,...
%     'EdgeColor',[0.95 0.95 0.95]);

%%
for x=1:size(ABK_B,2)
    ON(x) = size(find(ABK_B(:,x)>50),1) / reps * 100;
end

fig2 = figure('Name','% ON','NumberTitle','off');
set(fig2,'Position',[501 1 500 406]);                           hold on;
plot(12*ones(1,100),1:100,'--','Color',[1 0.69 0.39]);
plot(S_array,ON,'ob');
xlabel('N_S');                     ylabel('% ON');          
axis([S_array(1) S_array(end) 0 100]);                          hold off;
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');

%% Finish
close(pb);
clear progress pb a h i j r s t w x z tempA tempB leg* fig* p* temp;
clear Psa Psb Ta Tb Av Bv;
toc
