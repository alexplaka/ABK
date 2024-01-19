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

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear;      clc;     tic;                           
rng(1);

global agents k_b k_d k_s k_f k_r Km_f Km_r S;

totalTime = 1500;             % Simulation time (sec)
dt = 1/50;                   % Constant (fixed) time step increment (sec)
t_steps = totalTime / dt;
% time = zeros(1,t_steps);

agents = 100;
k_b = 0.02;                        % 0th order R synthesis rate, UPregulated by Ep
k_d = 0.075;                       % 1st order degradation rate (units: 1/sec)
k_s = 0.05;                        % 0th order R synthesis rate UPregulated by S
k_f = 1;                           % basal forward rate (1st order)
k_r = 0.05;                        % basal reverse rate (1st order)
Km_f = 5;                          % MICROSCOPIC Michaelis-Menten constant for forward rx
Km_r = 10;                         % MICROSCOPIC Michaelis-Menten constant for reverse rx

S = 15;                          % Assume number of S molecules/agents is NOT changing

% Initial population sizes
Ri = 0;        Epi = 0;      Ei = agents - Epi;					

[R_ss, Ep_ss] = RateCurve;           disp(['All R_ss values = ' num2str(R_ss')]);
% bifurc(Tr(1),Tep(1));

%% Symbolic calculations
% syms Ep R positive;
% 
% dR_sym  = k_b * Ep + k_s * S - k_d * R;
% dEp_sym = k_r * R * (agents - Ep) / (Km_r + agents - Ep) - k_f * Ep / (Km_f + Ep);
% 
% R_nc_sym = solve(dR_sym == 0,R,'MaxDegree',4,'Real',true);
% Ep_nc_sym = solve(dEp_sym == 0,Ep,'MaxDegree',4,'Real',true);
% 
% temp = subs(R_nc_sym,Ep,Ep_nc_sym(1));
% R_ss = double(solve(temp==R,R,'MaxDegree',4,'Real',true));
% % R_ss = double(vpasolve(temp==R,R,50));          % last argin is Initial Guess.
% Ep_ss = double(subs(Ep_nc_sym,R,R_ss));                      % clear temp;
% 
% Jac = jacobian([dR_sym,dEp_sym],[R,Ep]);              % Calculate Jacobian matrix
% 
% fprintf(1,'Steady-States:\t\t\tEigenvalues\n');
% for w=1:size(Ep_ss,1)
%     J = double(subs(Jac,[R,Ep],[R_ss(w),Ep_ss(w)]));
%     TrJ = J(1,1) + J(2,2);
%     DetJ = det(J);
%     DiscrJ = TrJ^2 - 4*DetJ;
%     lambda_plus(w) = (TrJ + sqrt(DiscrJ)) / 2;
%     lambda_minus(w) = (TrJ - sqrt(DiscrJ)) / 2;   
%     fprintf(1,['R=' num2str(R_ss(w)) ', Ep=' num2str(Ep_ss(w)) ':\t' ...
%         num2str(lambda_plus(w)) ', ' num2str(lambda_minus(w)) '\n']);
% end

%%
Tr = zeros(1,t_steps);
Tep = zeros(1,t_steps);             Te = zeros(1,t_steps);
P_f = zeros(1,t_steps);             P_r = zeros(1,t_steps);
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

%% ABK simulation
t = 1;

P_d = k_d * dt;                 % Probability of R degradation process (1st order wrt R)
P_s = k_s *S* dt;               % Probability of R synthesis process (0th order)

while t <= t_steps % && Tep(t)<agents
    P_b(t) = k_b * Tep(t) * dt;      % Probability of R synthesis process (0th order)
    P_f(t) = k_f / (Km_f + Tep(t)) * dt;          % wrt each Ep molecule
    P_r(t) = k_r * Tr(t) / (Km_r + Te(t)) * dt;   % wrt each E molecule
    
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
clear R Ep E;
finaltime = (t-1) * dt;
time = 0:dt:finaltime;
%% Solve DE
[t_sol, y_sol] = ode45(@mAct_dif,0:totalTime/500:totalTime,[Ri; Epi; Ei]);
%% Plot time course
figure('Name','Time Course','NumberTitle','off'); 
scatter(time,Tr,3,'ob');                                hold on;
scatter(time,Tep,3,'or');
% scatter(time,Te,3,'og');
scatter(t_sol,y_sol(:,1),3,'.b');                                       
scatter(t_sol,y_sol(:,2),3,'.r');
% scatter(t_sol,y_sol(:,3),3,'.g');                                     
axis([0 finaltime 0 agents]);                           hold off;
xlabel('time');                 ylabel('#');      
% legend('R stoc','Ep stoc','E stoc','R deter','Ep deter','E deter','Location','Best');

%% Plot Nullclines, State-space trajectories
syms R Ep positive;
dR_sym = + k_b * Ep + k_s * S - k_d * R;
dEp_sym = + k_r * (agents - Ep) * R / (Km_r + (agents - Ep)) - k_f * Ep / (Km_f + Ep);
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

[R_int Ep_int] = intersections(R_nc,array,array,Ep_nc(1,:));

[Rm,Epm] = meshgrid(0:agents/10:agents);     % Mesh-grid values for constructing state space
Em = agents - Epm;
dR  = + k_b .* Epm + k_s * S - k_d .* Rm;
dEp = + k_r .* Em .* Rm ./ (Km_r + Em) - k_f .* Epm ./ (Km_f + Epm);
r = ( dR.^2 + dEp.^2 ).^0.5;
% r=1;

figure('Name','State Space','NumberTitle','off');               hold on;
plot(R_nc,array,'b-.',array,Ep_nc,'m-.','LineWidth',1);
quiver(Rm,Epm,dR./r,dEp./r,'g');                  axis([0 agents 0 agents]);                  
scatter(y_sol(:,1),y_sol(:,2),'+k');       % plot deterministic trajectory
scatter(Tr,Tep,3,'.r');                    % plot stochastic trajectory
plot(R_int,Ep_int,'ob','MarkerSize',6);
xlabel('R');                    ylabel('E_p');                  hold off; 
legend('R nullcline','Ep nullcline','Direction field','Deterministic','ABK stochastic');
%% Finish
clear array R Ep dR dEp R_nc Ep_nc R_nc_sym Ep_nc_sym;
clear h i j k n r t temp tempR tempEp tempE;
toc

% - Works!!
% - Evidence that near a bifurcation point, stochasticity can lead to different switch states,
% in opposition to the deterministic prediction!!