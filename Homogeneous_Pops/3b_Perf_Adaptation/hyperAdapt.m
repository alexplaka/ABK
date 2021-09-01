% Hyperbolic Response Motif --> Adaptation
% 
%  S
% --> X -->              Rates: k_b*S, k_d
%     |
%     / 
% R <---> Rp             Rates: k_f, k_r
%     /
%     |
%     S                 S: Signal (reactant in forward reaction)

% Assume number of S molecules does not change.

% k_f: synthesis of Rp, 1st order process wrt R (upregulated by S), or 2nd order wrt R,S
% k_r: synthesis of R, 2nd order process (X, Rp), but X is not used up.
% k_b: synthesis of X, 0th order process (upregulated by S)
% k_d: degradation of X, 1st order process (X)

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;             tic;
clc; 
rng(0);
%% Declare variables and functions
global k_b k_d k_f k_r S;

maxTime = 40;                          % Maximum Simulation time (sec)
dt = 1/100;                            % Constant (fixed) time step increment (sec)
t_steps = maxTime / dt;

agents = 100;          	disp(['Agents = ' num2str(agents)]);
k_b = 0.1;                             % rate of synthesis of X; 0th order [units: 1/sec]
k_d = 0.1;                             % Degradation of X; 1st order [units: 1/sec]
k_f = 0.01;                            % basal forward rate (1st order)
k_r = 0.1;                             % basal reverse rate (2nd order)

S = 50;                                % Number of S molecules is constant

%% Simple formula for Rp_ss
phi = k_f * k_d / (k_r * k_b);
Rp_ss = phi / (1+phi) * agents;

%% Construct Rate curve

Rp_array = 0:agents;

dRpdt_degr = k_r .* Rp_array;
dRpdt_synt = k_f * S .* (agents - Rp_array);  

[Rp_ss_int, ssr] = intersections(Rp_array,dRpdt_degr,Rp_array,dRpdt_synt,0);
% Rp_int: Steady-state value of Rp derived from intersections in Rate curve
% ssr: Steady-state rate of forward and back reactions
figure('Name','Rate Curve','NumberTitle','off');
plot(Rp_array,dRpdt_degr,'r',Rp_array,dRpdt_synt,'-k',Rp_ss_int,ssr,'ob');
xlabel('R_p');        ylabel('dR_p/dt');    
legend('degradation','synthesis','R_{p-ss}','Location','NorthEast');
clear Rp_array;

%% Initialize - Preallocate memory for variables; ** Initial conditions **
t = 1;                              % Time counter variable
time = zeros(1,t_steps);            % For storing actual time values.

Tr = zeros(1,t_steps);              % Total sum of R agents in each time step	
Trp = zeros(1,t_steps);             % Total sum of Rp agents in each time step
Tx = zeros(1,t_steps);              % Total sum of R agents in each time step	

% ******** Initial conditions - Number of R, Rp, X Agents ********
Tr(1) = floor(agents);              Trp(1) = agents - Tr(1);
Tx(1) = 0;
% ************************************************************

tempR = zeros(1,agents);            tempX = zeros(1,agents);                    
% Put "1" where agenst are "alive", then randomize the array
for c=1:Tr(1),                      tempR(c)=1;             end
tempR = RandArray(tempR);           % Randomize array
Rv = [tempR ; tempR];               % Initialize vector storing state of R agents     
Rpv = ~Rv;                          % Initialize vector storing state of Rp agents 
for d=1:Tx(1),                      tempX(d)=1;             end
tempX = RandArray(tempX);
Xv = [tempX ; tempX];               % Initialize vector storing state of X agents
% Notes on Rv, Rpv:
% - Markov process, so only previous and current time steps needed --> 2 rows:          
% - Rv and Rpv are complementary (R+Rp=agents)
P_r = zeros(1,t_steps);
%% ABK Simulation

P_b = (k_b * S) * dt;               % P of X synthesis (0th order)

P_d = 1 - exp(-k_d * dt);           % P_ber of X degradation (1st order wrt X)
% P4 = k_d * dt;                     % P_dif of X degradation (1st order wrt X)

P_f = 1 - exp(-k_f * S * dt);    % P_ber for forward reaction: 1st order (wrt R)           
%     P_f = k_f * S * dt;              % P_dif for forward reaction: 1st order (wrt R)

% while abs(Trp(t)-Rp_ss) > 0.5 && t*dt <= maxTime
while t*dt <= maxTime
    
%     dt = 1 / (ko*Sa(t-1));            % Variable time step increment ; OLD
%     dt = exprnd(1/(k*S(t-1)));   % Exponentially-distributed time step increment ; OLD
    
%     P_r(t) = 1 - exp(-k_r * Tx(t) * dt);        % P_ber wrt Rp for reverse reaction: 2nd order            
    P_r(t) = k_r * Tx(t) * dt;                  % P_dif wrt Rp for reverse reaction: 2nd order

%   Take care of 0th order process first  
    if rand < P_b
        temp = find(Xv(1,:)==0);             % Randomly choose X agent which becomes "alive" 
        Xv(2,temp(ceil(rand * size(temp,2)))) = 1;  
    end
  
%   1st order degradation of X
    tempX = find(Xv(1,:)==1);
    for h = 1:size(tempX,2)                        % Reaction for each X molecule/agent
        if rand < P_d
            Xv(2,tempX(h)) = 0;                    % Degradation of X                         
        end
    end

    for i = 1:agents
        if Rv(1,i) == 1                     % if R agent is alive
            if rand < P_f                  % **check probability condition**
                Rv(2,i) = 0;                 % R agent "dies"
                Rpv(2,i) = 1;                 % R is converted to Rp
            end
        elseif Rpv(1,i) == 1					% if Rp agent is "alive"
            if rand < P_r(t)                  % **check probability condition**
                Rpv(2,i) = 0;                 % Rp agent "dies"
                Rv(2,i) = 1;                 % Rp is converted to R
            end
        end
    end
    
    Tr(t+1) = sum(Rv(2,:));      Trp(t+1) = sum(Rpv(2,:));       Tx(t+1) = sum(Xv(2,:));
    time(t+1) = time(t) + dt;
    
    Rv(1,:) = Rv(2,:);           Rpv(1,:) = Rpv(2,:);            Xv(1,:) = Xv(2,:);
    
    t = t + 1; 
end

disp(['ABK sim terminal R value   = ' num2str(Tr(end))]);  
disp(['ABK sim terminal Rp value  = ' num2str(Trp(end))]);  
disp(['ABK sim terminal X value   = ' num2str(Tx(end))]);  

% finaltime = (t-1) * dt;

%% Solve differential equation (was waiting for variable 'finaltime')
[t_sol, y_sol] = ode45(@hyperAdapt_dif,0:maxTime/100:maxTime,[Tr(1) ; Trp(1); Tx(1)]);

disp(['Diff eq terminal R  value   = ' num2str(y_sol(end,1))]);
disp(['Diff eq terminal Rp value   = ' num2str(y_sol(end,2))]);
disp(['Diff eq terminal X  value   = ' num2str(y_sol(end,3))]);

%% Graph ABK (stochastic) and deterministic results
fig1 = figure('Name','Time Course - Adaptation','NumberTitle','off');
set(fig1,'Position',[1 1 500 406]);                                   hold on;

% plot(time,Tr,'m','MarkerSize',1);              
plot(time,Trp,'b','MarkerSize',1);
plot(time,Tx,'Color',[0 0.8 0],'MarkerSize',1);

% plot(t_sol,y_sol(:,1),'--','Color',[1 0.75 0],'LineWidth',1);
plot(t_sol,y_sol(:,2),'--c','LineWidth',2);
plot(t_sol,y_sol(:,3),'--g','LineWidth',1);    

axis([0 maxTime 0 agents]);         
xlabel('t (sec)');        ylabel('N(t)');                              hold off;
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
% leg1 = legend('N_R(t)','N_{Rp}(t)','N_X(t)','DE N_R(t)',...
%     'DE N_{Rp}(t)','DE N_X(t)');                                         % With DE curves
leg1 = legend('N_{Rp}(t)','N_X(t)');                         % Without DE curves
set(leg1,'Location','NorthEast');
set(leg1,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95]);

%% Finish
clear c d h i j k r Rv Rpv Xv temp*;
toc
