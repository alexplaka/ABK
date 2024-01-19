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

% Considering individual agents.

% **************************************************************
% Considering heterogeneity in the regulation constant Km_r.
% **************************************************************

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear;      close all;   clc;     tic;                           
rng(1);

global agents k_b k_d k_s k_f k_r Km_f Km_r_mean S;

totalTime = 1500;                   % Simulation time (sec)
dt = 1/100;                         % Constant (fixed) time step increment (sec)
t_steps = totalTime / dt;
% time = zeros(1,t_steps);

agents = 100;
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

S = 15;                            % Assume number of S molecules/agents is NOT changing

% Initial population sizes
Ri = 0;        Epi = 0;      Ei = agents - Epi;					

[R_ss, Ep_ss] = RateCurve;           disp(['All R_ss values = ' num2str(R_ss')]);
% bifurc(Ri,Epi);

%% Set up simulation variables
Tr = zeros(1,t_steps+1);
Tep = zeros(1,t_steps+1);             Te = zeros(1,t_steps+1);

P_f = zeros(1,t_steps);               P_r = zeros(agents,t_steps);
P_b = zeros(1,t_steps);

avg_Km_r = zeros(1,t_steps+1);

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

tempE_log = E(1,:)==1;    % Logical array of existing E agents at time 0
avg_Km_r(1) = tempE_log * Km_r' / sum(tempE_log);    % Find initial <Km_r>

clear c d temp*;

%% ABK simulation
t = 1;

P_d = k_d * dt;                 % Probability of R degradation process (1st order wrt R)
P_s = k_s *S* dt;               % Probability of R synthesis process (0th order)

while t <= t_steps % && Tep(t)<agents
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
    
% The following implementation also works (treating it as a 1st order rx wrt S)     
%     for h=1:S                               % Reaction for each S molecule/agent
%         if rand < P_s
%             tempR = find(R(1,:)==0);         % Randomly choose R agent which becomes "alive" 
%             R(2,tempR(ceil(rand * size(tempR,2)))) = 1;
%         end
%     end
    
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
    
    tempE_log = E(2,:)==1;     % Logical array of existing E agents at end of time step
    avg_Km_r(t+1) = tempE_log * Km_r' / sum(tempE_log);    % Find <Km_r> at end of time step
    
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
% if t < t_steps
%     Tr = Tr(1:t);                                         
%     Tep = Tep(1:t);                             Te = Te(1:t);             
% end

clear R Ep E;

finaltime = (t-1) * dt;
time = 0:dt:finaltime;
%% Solve DE
[t_sol, y_sol] = ode45(@mAct_dif,0:totalTime/500:totalTime,[Ri; Epi; Ei]);

%% Plot time course
fig1 = figure('Name','Time Course','NumberTitle','off');               hold on;
set(fig1,'PaperPosition',[0 0 6 5],'PaperUnits','inches',...
    'WindowStyle','Normal','Position',[1 1 500 450]);

p(1) = plot(time(1:20:end),Tep(1:20:end),'b',...
    'DisplayName','$\textrm{ABK } N_{Ep}(t)$','LineWidth',0.5);
p(2) = plot(time(1:20:end),Tr(1:20:end),'r',...
    'DisplayName','$\textrm{ABK } N_R(t)$','LineWidth',0.5);   
p(3) = plot(t_sol,y_sol(:,2),'--c','DisplayName','$\textrm{DE } N_{Ep}(t)$','LineWidth',2);
p(4) = plot(t_sol,y_sol(:,1),'--m','DisplayName','$\textrm{DE } N_R(t)$','LineWidth',2); 

p(5) = plot(time,avg_Km_r,'k','DisplayName','$< \! K_{M,r}(t) \! >$','LineWidth',2);

axis([0 finaltime 0 agents]);                           hold off;
xlabel('$t \; \textrm{(sec)}$','Interpreter','LateX','FontSize',11);
ylabel('$N(t)$','Interpreter','LateX','FontSize',11);      
leg1 = legend(p);
set(leg1,'Location','NorthEast','FontSize',10,...
    'Interpreter','LateX','EdgeColor',[0.95 0.95 0.95]);

%% Plot time course FOR REPORT
fig1a = figure('Name','Time Course','NumberTitle','off');               hold on;
set(fig1a,'PaperPosition',[0 0 6 5],'PaperUnits','inches',...
    'WindowStyle','Normal','Position',[1 1 500 450]);

p(1) = plot(time(1:20:end),Tr(1:20:end),'r',...
    'DisplayName','$\textrm{ABK } N_R(t)$','LineWidth',0.5);   
p(2) = plot(time(1:20:end),Tep(1:20:end),'Color',[0 0.75 0],...
    'DisplayName','$\textrm{ABK } N_{Ep}(t)$','LineWidth',0.5);

p(3) = plot(time(1:20:end),Te(1:20:end),'Color',[0 0.5 0.75],...
    'DisplayName','$\textrm{ABK } N_{E}(t)$','LineWidth',0.5);

p(4) = plot(time,avg_Km_r,'k','DisplayName','$< \! K_{M,r}(t) \! >$','LineWidth',2);

% p_de(1) = plot(t_sol,y_sol(:,2),'--c','DisplayName','$\textrm{DE } N_{Ep}(t)$','LineWidth',2);
% p_de(2) = plot(t_sol,y_sol(:,1),'--m','DisplayName','$\textrm{DE } N_R(t)$','LineWidth',2); 

% Select one of the two title options below:
% title('$\textrm{Switch state OFF: Low } N_R \textrm{, low } N_{Ep} \textrm{, high } N_E$',...
%     'Interpreter','LateX','FontSize',12);
title('$\textrm{Switch state ON: High } N_R \textrm{, High } N_{Ep} \textrm{, low } N_E$',...
    'Interpreter','LateX','FontSize',12);

axis([0 finaltime 0 agents]);                           hold off;
xlabel('$t \; \textrm{(sec)}$','Interpreter','LateX','FontSize',11);
ylabel('$N(t)$','Interpreter','LateX','FontSize',11);     

leg1a = legend(p);
set(leg1a,'Location','NorthEast','FontSize',10,...
    'Interpreter','LateX','EdgeColor',[0.95 0.95 0.95]);


%% Compute Nullclines
syms R Ep positive;
dR_sym = + k_b * Ep + k_s * S - k_d * R;
dEp_sym = + k_r * (agents - Ep) * R / (Km_r_mean + (agents - Ep)) - k_f * Ep / (Km_f + Ep);
% Nullclines
R_nc_sym = solve(dR_sym == 0,R);
Ep_nc_sym = solve(dEp_sym == 0,Ep);

array = 0:agents;    
if k_f/k_r < agents
    array(k_f/k_r + 1) = [];      % Remove R=k_f/k_r entry (avoid div by 0)
end
R_nc  = double(subs(R_nc_sym,Ep,array));
Ep_nc = double(subs(Ep_nc_sym,R,array));

if size(R_nc,1)>1                % Nullcline values must be positive!
    for j=size(R_nc,1):-1:1
        if isempty(find(sign(R_nc(j,:)) == -1, 1,'first')) == false
            R_nc(j,:) = [];
        end
    end
end

if size(Ep_nc,1)>1               % Nullcline values must be positive!
    for k=size(Ep_nc,1):-1:1
        if isempty(find(sign(Ep_nc(k,:)) == -1, 1,'first')) == false
            Ep_nc(k,:) = [];
        end
    end
end

[R_int, Ep_int] = intersections(R_nc,array,array,Ep_nc(1,:));

[Epm, Rm] = meshgrid(0:agents/10:agents);     % Mesh-grid values for constructing state space
Em = agents - Epm;
dR  = + k_b .* Epm + k_s * S - k_d .* Rm;
dEp = + k_r .* Em .* Rm ./ (Km_r_mean + Em) - k_f .* Epm ./ (Km_f + Epm);
% r = ( dR.^2 + dEp.^2 ).^0.5;
r=1;

%% Now plot phase plane
fig2 = figure('Name','State Space','NumberTitle','off');               hold on;
set(fig2,'PaperPosition',[0 0 6 5],'PaperUnits','inches',...
    'WindowStyle','Normal','Position',[1 1 500 450]);

p2(1) = plot(array,R_nc,'LineWidth',2,'Color',[0.4 0.4 1],...
    'DisplayName','$N_R \textrm{ nullcline}$');
p2(2) = plot(Ep_nc,array,'LineWidth',2,'Color',[1 0.3 0.2],...
    'DisplayName','$N_{Ep} \textrm{ nullcline}$');

% p2(3) = quiver(Epm,Rm,dEp./r,dR./r,'g');                     
p2(3) = streakarrow(Epm,Rm,dEp./r,dR./r,0.7,1);             % Direction field
set(p2(3),'DisplayName','Direction Field');

% plot deterministic trajectory
p2(4) = plot(y_sol(:,2),y_sol(:,1),'+k','DisplayName','DE  Trajectory','MarkerSize',8);

plot(Tep,Tr,'.r');                    % plot stochastic trajectory

if size(R_int,1)==3
    for z=1:3
        if mod(z,2)==0              % Unstable fixed point
            p_fp(z) = plot(Ep_int(z),R_int(z),'oc','MarkerSize',8,...
                'MarkerFaceColor','r','DisplayName','Saddle FP');        
        else                        % Two Stable fixed points
            p_fp(z) = plot(Ep_int(z),R_int(z),'oc','MarkerSize',8,...
                'MarkerFaceColor','b','DisplayName','Stable FP');       
        end
    end
else
    p_fp = plot(Ep_int,R_int,'oc','MarkerSize',8,...
        'MarkerFaceColor','b','DisplayName','Stable FP');          % Sole stable fixed point
end

axis([0 agents 0 agents]);   
ylabel('$N_R$','Interpreter','LateX','FontSize',11);                    
xlabel('$N_{Ep}$','Interpreter','LateX','FontSize',11);                  hold off; 

leg2 = legend([p2(2:-1:1), p_fp(1:2)]);
set(leg2,'Location','NorthWest','FontSize',10,...
    'Interpreter','LateX','EdgeColor',[0.95 0.95 0.95]);

%% Finish
clear array R Ep dR dEp R_nc Ep_nc R_nc_sym Ep_nc_sym;
clear h i j k n r t z temp* p* leg* fig*;
toc
