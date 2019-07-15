%   Ep  <--->  E            Rates: k_f, k_r (Reverse rx activated by R)
%    |   ^                  Michaelis-Menten constants: Km_f, Km_r
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
% k_f: E synthesis, MM process
% k_r: Ep synthesis, MM process wrt R (but R is not consummed)
% Km_f, Km_r: MM constants for forward and reverse rxs

% Note: Assume number of S molecules/agents is NOT changing (eg, S is DNA molecules)

% Here, we construct the SR curve for this motif.

% **************************************************************
% Considering heterogeneity in the regulation constant Km_r.
% **************************************************************

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;      close all;       clc;     tic;                           % rng(1);

global agents k_b k_d k_s k_f k_r Km_f Km_r_mean S;

totalTime = 1500;             % Simulation time (sec)
dt = 1/100;                   % Constant (fixed) time step increment (sec)
t_steps = totalTime / dt;
% time = zeros(1,t_steps);

agents = 100;                      % Max number of agents for any species (choose even number)
k_b = 0.02;                        % 0th order R synthesis rate, UPregulated by Ep
k_d = 0.075;                       % 1st order degradation rate (units: 1/sec)
k_s = 0.05;                        % 0th order R synthesis rate UPregulated by S
k_f = 1;                           % basal forward rate (1st order)
k_r = 0.05;                        % basal reverse rate (1st order)
Km_f = 5;                          % MICROSCOPIC Michaelis-Menten constant for forward rx
Km_r_mean = 10;                    % MICROSCOPIC Michaelis-Menten constant for reverse rx

% Agent-specific Km_r values for species E. Comment one of the following two options:
% Option 1: Random integers in 1:19 (mean=10, sdev=5.5)
% Km_r = randi(19,1,agents);  
% Option 2: Two subpopulations of equal size with Km_r=1,19 (mean=10, sdev=9)
Km_r = [ones(1,agents/2) 19*ones(1,agents/2)];

S_array = 0:30;

% Initial population sizes
Ri = 0;        Epi = 0;      Ei = agents - Epi;					

% Preallocate memory for ABK_R_ss, and time series data
ABK_R_ss = zeros(1, size(S_array,2));
Rv = zeros(size(S_array,2), t_steps+1);
Epv = zeros(size(S_array,2), t_steps+1);


for s = 1:size(S_array,2)
    
    fprintf(1,'.');
    
    S = S_array(s);                    % Assume number of S molecules/agents is NOT changing
    
    Tr = zeros(1,t_steps);
    Tep = zeros(1,t_steps);             Te = zeros(1,t_steps);
    P_f = zeros(1,t_steps);             P_r = zeros(agents,t_steps);
    P_b = zeros(1,t_steps);

    % ******** Initial conditions - Number of R, Ep, E Agents ********
    Tr(1) = Ri;             Tep(1) = Epi;           Te(1) = Ei;					
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

    while t <= t_steps 
        P_b(t) = k_b * Tep(t) * dt;      % Probability of R synthesis process (0th order)
        P_f(t) = k_f / (Km_f + Tep(t)) * dt;          % wrt each Ep molecule

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
            P_r(tempE(j),t) = k_r * Tr(t) / ( Km_r(tempE(j)) + Te(t) ) * dt;   % wrt each E molecule
            if rand < P_r(tempE(j),t)
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

    Rv(s,:) = Tr;                           Epv(s,:) = Tep;

    %   Average R population over the last 750 sec of the simulation
    %   Steady-state has been reached for this time period
    ABK_R_ss(s) = mean(Tr(end-750/dt:end));

end                 % end 'for s'

finaltime = (t-1) * dt;
time = 0:dt:finaltime;

%% Plot signal-response (SR) curve
fig1 = figure('Name','Signal-Response Curve','NumberTitle','off');
set(fig1,'Position',[1 1 500 406]);                                 hold on;

% Plot Theoretical SR curve
SR = load('./Theory/SR_curve.mat','S_array','R_ss');
for h=1:size(SR.R_ss,2)
    if mod(h,2) == 0
        p(h) = plot(SR.S_array,SR.R_ss(:,h),'--','LineWidth',2,...
        'Color',[1 0.69 0.39],'DisplayName','Unstable');
    else
    p(h) = plot(SR.S_array,SR.R_ss(:,h),'-b','LineWidth',2,'DisplayName','Stable');
    end
end

% Plot ABK results
p_ABK = plot(S_array,ABK_R_ss,'.r','DisplayName','ABK');
xlabel('N_S');                               
ylabel('N_R^*');                            

axis([0 S_array(end) 0 max(ABK_R_ss)+5]);                                    hold off;
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
leg1 = legend([p(1), p(2), p_ABK],'Location','SouthEast');
set(leg1,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95]);

%% Finish - Notes
clear h i j k n t Tr Tep Te R Ep E temp tempR tempEp tempE;
toc
