%   E   <--->  Ep             Rates: k_r, k_f 
%    |     ^                  Michaelis-Menten constants: Km_f, Km_r
%    |    /                   Reverse rx E --> Ep activated by R
%    \    |   k_d
%    -->  R   -->                      R: Response
%   k_b        ^   
%              |
%              S                       S: Signal

% Simulating homeostatic process for R: 
% k_b: R synthesis, 1st order process wrt Ep (but Ep is not consummed)
% k_d: R degradation, 2nd order process wrt R, S
% k_f: E synthesis, MM process
% k_r: Ep synthesis, MM process wrt R (but R is not consummed)
% Km_f, Km_r: MM constants for forward and reverse rxs

% Note: Assume number of S molecules/agents is NOT changing (eg, S is DNA molecules)
% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;      clc;     tic;                           % rng(1);
rng(0);

global maxRagents agents k_b k_d k_f k_r Km_f Km_r S;

totalTime = 1000;             % Simulation time (sec)
dt = 1/50;                   % Constant (fixed) time step increment (sec)
t_steps = totalTime / dt;
% time = zeros(1,t_steps);

agents = 100;                      % Max # of E, Ep agents
maxRagents = 100;                  % Max # of R agents
k_b = 0.02;                        % 1st order R synthesis rate wrt E
k_d = 0.002;                       % MICROSCOPIC 2nd order degradation rate (units: 1/sec)
k_f = 1;                           % basal forward rate (Ep --> E)
k_r = 0.05;                        % basal reverse rate (E --> Ep)
Km_f = 5;                          % MICROSCOPIC Michaelis-Menten constant for forward rx
Km_r = 10;                         % MICROSCOPIC Michaelis-Menten constant for reverse rx
S = 10;                            % Assume number of S molecules/agents is NOT changing

Tr = zeros(1,t_steps);
Tep = zeros(1,t_steps);             Te = zeros(1,t_steps);
P_f = zeros(1,t_steps);             P_r = zeros(1,t_steps);
P_b = zeros(1,t_steps);

% ******** Initial conditions - Number of R, Ep, E Agents ********
Tr(1) = 20;        Tep(1) = 30;      Te(1) = agents - Tep(1);					
% ****************************************************************

tempR = zeros(1,maxRagents);        tempEp = zeros(1,agents);             
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

[R_ss, E_ss] = RateCurve;           disp(['All R_ss values = ' num2str(R_ss')]);
%% ABK simulation
t = 1;
P_d = k_d * S * dt;          % Probability of R degradation wrt R

while t < t_steps && Tr(t) < maxRagents
    P_b(t) = k_b * Te(t) * dt;   % Probability of R synthesis: 0th order UPreg by E
    P_f(t) = k_f / (Km_f + Tep(t)) * dt;          % wrt each Ep molecule
    P_r(t) = k_r * Tr(t) / (Km_r + Te(t)) * dt;   % wrt each E molecule

%   Take care of 0th order processes first; 
%   treat R synthesis as 0th order, UPregulated by E (b/c E is not consumed)
    if rand < P_b(t)
        tempR = find(R(1,:)==0);    % Randomly choose R agent synthesis 
        R(2,tempR(ceil(rand * size(tempR,2)))) = 1;
    end
%   End of 0th order processes
    
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
    P_b = P_b(1:t);
end
clear R Ep E;
finaltime = (t-1) * dt;
time = 0:dt:finaltime;
% Solve DE
[t_sol, y_sol] = ode45(@homeo_dif,0:finaltime/500:finaltime,[Tr(1); Tep(1); Te(1)]);
%% Plot time course
figure('Name','Time Course','NumberTitle','off');           hold on;
scatter(time,Tr,3,'ob');                                
scatter(time,Tep,3,'or');
scatter(time,Te,3,'og');
scatter(t_sol,y_sol(:,1),3,'.b');                                       
scatter(t_sol,y_sol(:,2),3,'.r');
scatter(t_sol,y_sol(:,3),3,'.g');                                     
axis([0 finaltime 0 maxRagents]);                           hold off;
xlabel('time');                 ylabel('#');      
legend('R stoc','Ep stoc','E stoc','R deter','Ep deter','E deter','Location','Best');

%% Plot Nullclines, State-space trajectories
syms R E positive;
dR_sym = + k_b * E - k_d * R * S;
dE_sym = - k_r * E * R / (Km_r + E) + k_f * (agents - E) / (Km_f + (agents - E));
% Nullclines
R_nc_sym = solve(dR_sym == 0,R);
E_nc_sym = solve(dE_sym == 0,E);
array = 0:maxRagents;           
array(k_f/k_r + 1) = [];      % Remove R=k_f/k_r entry (avoid div by 0)
R_nc  = double(subs(R_nc_sym,E,array));
E_nc = double(subs(E_nc_sym,R,array));

if size(R_nc,1)>1                % Nullcline values must be positive!
    for j=size(R_nc,1):-1:1
        if isempty(find(sign(R_nc(j,:)) == -1, 1,'first')) == false
            R_nc(j,:) = [];
        end
    end
end

if size(E_nc,1)>1               % Nullcline values must be positive!
    for k=size(E_nc,1):-1:1
        if isempty(find(sign(E_nc(k,:)) == -1, 1,'first')) == false
            E_nc(k,:) = [];
        end
    end
end

[R_int, E_int] = intersections(R_nc,array,array,E_nc);

% Mesh-grid values for constructing state space
[Rm,Em] = meshgrid(0:maxRagents/10:maxRagents,0:agents/10:agents);  
dR  = + k_b .* Em - k_d .* Rm * S;
dE = - k_r .* Em .* Rm ./ (Km_r + Em) + k_f .* (agents-Em) ./ (Km_f + (agents-Em));
r = ( dR.^2 + dE.^2 ).^0.5;
% r=1;

figure('Name','State Space','NumberTitle','off');               hold on;
plot(R_nc,array,'b-.',array,E_nc,'m-.','LineWidth',1);
quiver(Rm,Em,dR./r,dE./r,'g');                  axis([0 maxRagents 0 agents]);                  
scatter(y_sol(:,1),y_sol(:,3),'+k');       % plot deterministic trajectory
scatter(Tr,Te,3,'.r');                     % plot stochastic trajectory
plot(R_int,E_int,'ob','MarkerSize',6);
xlabel('R');                    ylabel('E');                  hold off; 
legend('R nullcline','E nullcline','Direction field','Deterministic','ABK stochastic');

clear array R Ep dR dEp R_nc Ep_nc R_nc_sym Ep_nc_sym;
clear h i j k n r t temp tempR tempEp tempE;
toc
