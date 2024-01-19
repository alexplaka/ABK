%    R  <--->  Rp            Rates: k_f, k_r (Forward rx activated by X)
%          ^    |            Michaelis-Menten constants: Km_f, Km_r
%         /     |
%         |     \
%    -->  X    -->                      Rp: Response
%   k_b   ^    k_d1/k_d2
%         | k_s
%         S                             S: Signal

% Simulating negative feedback process: 
% k_b: X synthesis, 0th order process
% k_s: X synthesis, 1st order process wrt S (but S is not consummed)
% k_s: ALTERNATIVE IMPLEMENTATION: X synthesis, 0th order process, UPregulated by S
% k_d1: X degradation, 1st order process wrt X
% k_d2: X degradation, 2nd order process wrt X, Rp (but Rp is not consummed)
% k_f: Rp synthesis, MM process wrt X (but X is not consummed)
% k_r: R synthesis, MM process 
% Km_f, Km_r: MM constants for forward and reverse rxs

% Note: Assume number of S molecules/agents is NOT changing (eg, S is DNA molecules)

% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear;      clc;     tic;                           % rng(1);
rng(0);
global agents k_s k_b k_d1 k_d2 k_f k_r Km_f Km_r S;

totalTime = 1000;             % Simulation time (sec)
dt = 1/100;                   % Constant (fixed) time step increment (sec)
t_steps = totalTime / dt;
% time = zeros(1,t_steps);

agents = 100;
k_b = 0;                        % 0th order X synthesis rate
k_d1 = 0;                       % 1st order degradation rate (units: 1/sec)
k_d2 = 0.005;                       % 2nd order degradation rate wrt X, Rp (units: 1/sec)
k_s = 0.10;                         % 0th order R synthesis rate UPregulated by S
k_f = 0.01;                           % basal forward rate (1st order)
k_r = 1;                        % basal reverse rate (1st order)
Km_f = 10;                          % MICROSCOPIC Michaelis-Menten constant for forward rx
Km_r = 10;                           % MICROSCOPIC Michaelis-Menten constant for reverse rx

S = 30;                            % Assume number of S molecules/agents is NOT changing

Xi = 0;                             Rpi = 0;
%% Symbolic calculations - Linearize and calculate eigenvalues
syms X Rp positive;
dX_sym = + k_b + k_s * S - k_d1 * X - k_d2 * X * Rp;
dRp_sym = + k_f * (agents - Rp) * X / (Km_f + (agents-Rp)) - k_r * Rp / (Km_r + Rp);
% Nullclines
X_nc_sym = solve(dX_sym == 0,X,'MaxDegree',4,'Real',true);
Rp_nc_sym = solve(dRp_sym == 0,Rp,'MaxDegree',4,'Real',true);

%% Initialize, set Initial conditions.
Tx = zeros(1,t_steps);
Tr = zeros(1,t_steps);              Trp = zeros(1,t_steps);
P_f = zeros(1,t_steps);             P_r = zeros(1,t_steps);
P_d2 = zeros(1,t_steps);

% ******** Initial conditions - Number of R, Ep, E Agents ********
Tx(1) = Xi;         Trp(1) = Rpi;         Tr(1) = agents - Trp(1);					
% ****************************************************************

tempX = zeros(1,2*agents);        tempR = zeros(1,agents);             
% Put "1" where agents are "alive", then randomize the array
for c=1:Tx(1),                  tempX(c)=1;             end
for d=1:Tr(1),                  tempR(d)=1;            end
tempX = RandArray(tempX);                   % Randomize R array
tempR = RandArray(tempR);                   % Randomize R array
% Markov process, so only previous and current time steps needed --> 2 rows:    
Xv = [tempX ; tempX];                        % Initialize Vector storing state of X agents
Rv = [tempR ; tempR];                        % Initialize Vector storing state of R agents
Rpv = ~ Rv;                                  % Initialize Vector storing state of Rp agents
% - R and Rp are complementary (R + Rp = agents)
clear c d tempX tempR;
%% ABK simulation
t = 1;

P_b = k_b * dt;                  % Probability of X synthesis process (0th order)
P_s = k_s * S * dt;              % Probability of X synthesis process (0th order), UPreg'd by S
P_d1 = 1-exp(-k_d1 * dt);        % Probability of X degradation process (1st order) wrt X

while t <= t_steps % && Tep(t)<agents
    P_d2(t) = k_d2 * Trp(t) * dt;                          % wrt each X molecule
    P_f(t) = k_f * Tx(t) / (Km_f + Tr(t)) * dt;            % wrt each R molecule
    P_r(t) = k_r / (Km_r + Trp(t)) * dt;                   % wrt each Rp molecule
    
    % Take care of 0th order processes first
    if rand < P_b
        tempX = find(Xv(1,:)==0);                % Randomly choose X agent synthesis 
        Xv(2,tempX(ceil(rand * size(tempX,2)))) = 1;
    end
    
    if rand < P_s
        tempX = find(Xv(1,:)==0);                % Randomly choose X agent synthesis 
        Xv(2,tempX(ceil(rand * size(tempX,2)))) = 1;
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
            Rv(2,tempR(j)) = 0;         % Conversion of R ...
            Rpv(2,tempR(j)) = 1;        % to Rp
        end
    end
    
    tempRp = find(Rpv(1,:)==1);
    for k=1:size(tempRp,2)
        if rand < P_r(t)
            Rpv(2,tempRp(k)) = 0;       % Conversion of Rp ...
            Rv(2,tempRp(k)) = 1;        % to R
        end
    end
    
    Tx(t+1) = sum(Xv(2,:));
    Tr(t+1) = sum(Rv(2,:));             Trp(t+1) = sum(Rpv(2,:));      
 
    Xv(1,:) = Xv(2,:);
    Rv(1,:) = Rv(2,:);                  Rpv(1,:) = Rpv(2,:);
    t = t + 1;
end
% Remove unnecessary terminal 0's from arrays
if t < t_steps
    Tx = Tx(1:t);               Tr = Tr(1:t);           Trp = Trp(1:t);                                           
    P_f = P_f(1:t);             P_r = P_r(1:t);         P_d2 = P_d2(1:t);                 
end

%% Solve DE
% finaltime = (t-1) * dt;
finaltime = totalTime;
[t_sol, y_sol] = ode23(@negFb_dif,0:finaltime/500:finaltime,[Xi; Rpi]);
%% Plot time course
time = 0:dt:finaltime;
fig0 = figure('Name','Time Course','NumberTitle','off');           hold on;
set(fig0,'Position',[1 1 500 450]);
fig0.PaperUnits = 'inches';
fig0.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.

plot(time,Tx,'b');                                
plot(time,Trp,'r');
plot(t_sol,y_sol(:,1),'--','Color',[0.15 1 0.75],'LineWidth',2); 
plot(t_sol,y_sol(:,2),'--','Color',[0.75 0 1],'LineWidth',2);  
axis([0 finaltime 0 agents]);                               hold off;
xlabel('t (sec)');                 ylabel('N(t)');  

leg0 = legend('ABK N_R(t)','ABK N_{Xp}(t)','DE N_R(t)','DE N_{Xp}(t)','Location','Best');
set(leg0,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');

%% Calculate Nullclines, Direction Field
array = 1:agents;
k = k_r/k_f;
if k <= agents
    if mod(k,1) == 0 
        array(k) = [];      % Remove X=k_r/k_f entry (avoid div by 0 when calc Rp_nc)
    end
end
X_nc  = double(subs(X_nc_sym,Rp,array));
Rp_nc = double(subs(Rp_nc_sym,X,array));   % result has 2 entries, +/-

if size(X_nc,1)>1                % Nullcline values must be positive!
    for j=size(X_nc,1):-1:1
        if isempty(find(sign(X_nc(j,:)) == -1, 1,'first')) == false
            X_nc(j,:) = [];
        end
    end
end

if size(Rp_nc,1)>1               % Nullcline values must be positive!
    for k=size(Rp_nc,1):-1:1
        if isempty(find(sign(Rp_nc(k,:)) == -1, 1,'first')) == false
            Rp_nc(k,:) = [];
        end
    end
end
% If above does not eliminate one row or Rp_nc (i.e., all values are positive)...
Rp_nc = Rp_nc(2,:);

[X_int Rp_int] = intersections(X_nc,array,array,Rp_nc);

[Xm,Rpm] = meshgrid(0:5:agents);     % Mesh-grid values for constructing state space
Rm = agents - Rpm;
dX  = + k_b + k_s * S - k_d1 .* Xm - k_d2 .* Xm .* Rpm;
dRp = + k_f .* Rm .* Xm ./ (Km_f + Rm) - k_r .* Rpm ./ (Km_r + Rpm);
% r = ( dX.^2 + dRp.^2 ).^0.5;
r=1;

%% Plot Nullclines, State-space trajectories
% **** Using Report nomenclature X |-> R , Rp |-> Xp ****

fig1 = figure('Name','State Space','NumberTitle','off');               hold on;
set(fig1,'Position',[1 1 500 450]);                                                    
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.

% ** Plot nullclines **
p_ncR = plot(X_nc,array,'','Color',[1 0.3 0.2],'LineWidth',2,'DisplayName','N_R nullcline');
p_ncXp = plot(array,Rp_nc,'','Color',[0.4 0.4 1],'LineWidth',2,'DisplayName','N_{Xp} nullcline');      

% ** Plot Direction field ** : Choose 1 of the below options
% quiver(Xm,Rpm,dX./r,dRp./r,'g');                                   
streakarrow(Xm,Rpm,dX./r,dRp./r,0.7,1);             

% plot deterministic trajectory
p_de = plot(y_sol(:,1),y_sol(:,2),'+k','MarkerSize',3,'DisplayName','DE Trajectory');      

% plot full or partial stochastic trajectory
% p_abk = scatter(Tx,Trp,3,'.r');                         
% p_abk = scatter(Tx(1,250/dt:350/dt),Trp(1,250/dt:350/dt),3,'.r'); 

% ** Fixed pt ** : Choose 1 of the below options
p_fp = plot(X_int,Rp_int,'oc','MarkerSize',6,'MarkerFaceColor','b','DisplayName','Stable FP'); 
    
% axis([0 agents 0 agents]);  
axis([0 100 0 30]);                     % Custom axes ranges
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');               

% comet(Tx(1,250/dt:350/dt),Trp(1,250/dt:350/dt))

xlabel('R');                    ylabel('X_p');                  hold off;     

% legend('X nullcline','Rp nullcline','Direction field','DE','ABK','Stable FP');
leg1 = legend([p_ncR p_ncXp p_de p_fp]);
set(leg1,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');

%% Finish
clear array X Rp dX dRp X_nc Rp_nc X_nc_sym Rp_nc_sym;
clear h i j k n r t temp tempR tempRp tempX Xv Rv Rpv;
toc
