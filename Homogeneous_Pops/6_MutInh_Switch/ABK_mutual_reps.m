%    k_da           k_db
%    <--  A _  _ B  -->                 % For Report: B = R
%         ^  \/  ^
%     k_a |  /\  | k_b
%         | /  \ |
%         |~    ~|
% A and B are synthesized (0th order constants k_a, k_b respectively)
% A and B are degraded (k_da, k_db)
% A and B inhibit each other's synthesis (MUTUAL inhibition, alpha < 1):
% B influences the rate of the synthesis of A
% A influences the rate of the synthesis of B

% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;            tic;                  clc; 
rng('shuffle');
pb = waitbar(0,'0');
%% Declare variables and functions
global k_da k_db agents;
global k_a alpha_b K_b n_b;
global k_b alpha_a K_a n_a;

agents = 200;                       % disp(['Agents = ' num2str(agents)]);

maxTime = 2000;                        % Maximum Simulation time (sec)
dt = 1/50;                              % Constant (fixed) time step increment (sec)
t_steps = maxTime / dt;

reps = 10;                             % Number of times to repeat experiment

A_all = zeros(reps,t_steps+1);        B_all = zeros(reps,t_steps+1);

% ** PARAMETERS **
k_a = 1;                            % MICROSCOPIC basal rate of A synthesis (0th order)
k_b = 1;                            % MICROSCOPIC basal rate of B synthesis (0th order)
k_da = 0.01;                        % degradation rate constant (1st order) for A
k_db = 0.01;                        % degradation rate constant (1st order) for B
% Populations producing half-maximal regulatory effect
K_a = 50;                       K_b =  50;

% *** Feedback parameters ***
% Degree of activation: >1 activator, <1 repressor, =1 no regulation
alpha_a = 0;                        % Effect of A on synthesis of B; = 0 means Complete Repression
alpha_b = 0;                        % Effect of B on synthesis of A; = 0 means Complete Repression
% Hill Coefficient: measure of cooperativity
n_b = 3;                      n_a = 3;

% ******** Initial conditions - Number of A, B Agents ********
Ai = 30;                      Bi = 75;
% ************************************************************

%% Run repetitions of ABK simulation
for n = 1:reps
    
    progress = n / reps;
    waitbar(progress,pb,sprintf('%.0f%% %',progress*100));
    
    % Initialize - Preallocate memory for variables; ** Initial conditions **
    t = 1;                              % Time counter variable
    P_sa = zeros(1,t_steps);            % For storing probability of A synthesis at each time step
    P_sb = zeros(1,t_steps);            % For storing probability of B synthesis at each time step

    Ta = zeros(1,t_steps);              % Sum of A agents in each time step	
    Tb = zeros(1,t_steps);              % Sum of B agents in each time step

    % ******** Initial conditions - Number of A, B Agents ********
    Ta(1) = Ai;                         Tb(1) = Bi;					
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
%     clear c d tempA tempB;

    % ABK Simulation
    P_da = 1 - exp(-k_da * dt);                % P_ber of degradation of A, 1st order wrt A
%     P_da = k_da * dt;                     % P_dif of degradation of A, 1st order wrt A
    P_db = 1 - exp(-k_db * dt);                % P_ber of degradation of B, 1st order wrt B
%     P_db = k_db * dt;                     % P_dif of degradation of B, 1st order wrt B

    while t <= t_steps % && Ta(t) > agents/100
                
        k_sa = k_a * (K_b^n_b + alpha_b * Tb(t)^n_b) / (K_b^n_b + Tb(t)^n_b);
        k_sb = k_b * (K_a^n_a + alpha_a * Ta(t)^n_a) / (K_a^n_a + Ta(t)^n_a);
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
%     if t < t_steps
%         Ta = Ta(1:t);                            Tb = Tb(1:t);          
%     end
    % disp(['ABK sim terminal A value   = ' num2str(Ta(end))]);
    % disp(['ABK sim terminal B value   = ' num2str(Tb(end))]);                   clear Av Bv;
    
    % Save each trial's time trajectory
    A_all(n,:) = Ta;                           B_all(n,:) = Tb;  

end         % end 'for n' loop

%% Solve differential equation (was waiting for variable 'finaltime')
[t_sol, y_sol] = ode45(@mutual_dif,0:maxTime/100:maxTime,[Ai ; Bi]);
fprintf(1,['\nDiff eq terminal A value =\t' num2str(y_sol(end,1)) '\n']);
fprintf(1,['Diff eq terminal B value =\t' num2str(y_sol(end,2)) '\n']);

%% Calculate nullclines, direction field
% *** Set up differential equations of system using symbolic variables ***
syms A B k_sa_sym k_sb_sym positive;

k_sa_sym = k_a * (K_b^n_b + alpha_b * B^n_b) / (K_b^n_b + B^n_b);
k_sb_sym = k_b * (K_a^n_a + alpha_a * A^n_a) / (K_a^n_a + A^n_a);

dA_sym = + k_sa_sym - k_da * A;
dB_sym = + k_sb_sym - k_db * B;

% Nullclines
A_nc_sym = solve(dA_sym == 0,A,'MaxDegree',4,'Real',true);
B_nc_sym = solve(dB_sym == 0,B,'MaxDegree',4,'Real',true);

A_nc = double(subs(A_nc_sym,B,0:agents));
B_nc = double(subs(B_nc_sym,A,0:agents));

if size(B_nc,1) > 1                % Nullcline values must be positive!
    for j=size(B_nc,1):-1:1
        if isempty(find(sign(B_nc(j,:)) == -1, 1,'first')) == false
            B_nc(j,:) = [];
        end
    end
end

[Am,Bm] = meshgrid(1:agents/15:agents);     % Mesh-grid values for constructing state space
k_sa_num = subs(k_sa_sym,[A,B],{Am,Bm});
k_sb_num = subs(k_sb_sym,[A,B],{Am,Bm});
dA_num = double(subs(dA_sym,[k_sa_sym,k_sb_sym,A,B],{k_sa_num,k_sb_num,Am,Bm}));
dB_num = double(subs(dB_sym,[k_sa_sym,k_sb_sym,A,B],{k_sa_num,k_sb_num,Am,Bm}));
% r = ( dA_num.^2 + dB_num.^2 ).^0.5;
r = 1;
Ua = dA_num./r;          Ub = dB_num./r;
V = sqrt(Ua.^2 + Ub.^2);                % Magnitude of velocity vector

[A_ss_int, B_ss_int] = intersections(A_nc,0:agents,0:agents,B_nc);

% Clean up duplicate ss data
for z=size(B_ss_int,1):-1:2
    temp = B_ss_int(z) - B_ss_int(z-1);
    if abs(temp) < 0.0001
        B_ss_int(z) = [];
        A_ss_int(z) = [];
    end
end

Jac = jacobian([dA_sym,dB_sym],[A,B]);                  % Calculate Jacobian matrix

fprintf(1,'\nFixed Points:\t\tEigenvalues\n');
for w=1:size(A_ss_int,1)
    J = double(subs(Jac,[A,B],[A_ss_int(w),B_ss_int(w)]));
    TrJ = J(1,1) + J(2,2);
    DetJ = det(J);
    DiscrJ = TrJ^2 - 4*DetJ;
    lambda_plus(w) = (TrJ + sqrt(DiscrJ)) / 2;
    lambda_minus(w) = (TrJ - sqrt(DiscrJ)) / 2;   
    fprintf(1,['A= %4.2f , B= %4.2f :\t%+6.4f , %+6.4f \n'],...
        A_ss_int(w),B_ss_int(w),lambda_plus(w),lambda_minus(w));
end

%% Classify trajectories in terms of which equilibrium (fixed) point they reach
% Will average the population size for the last 1000 sec of each simulation (SS is reached by then).

sink1 = [];             sink2 = [];             % source = []; 

for n=1:reps
    
    term_avgA(n) = mean(A_all(n,1000/dt:end));          % Terminal avg population size of A
    term_avgB(n) = mean(B_all(n,1000/dt:end));          % Terminal avg population size of B

    if term_avgA(n) < 2.5*A_ss_int(1)
        sink1 = [sink1 n];
        
        A_all_sink1(size(sink1,2),:) = A_all(n,:);                   
        B_all_sink1(size(sink1,2),:) = B_all(n,:);                   
        
    elseif term_avgB(n) < 3.5*B_ss_int(3)
        sink2 = [sink2 n];
        
        A_all_sink2(size(sink2,2),:) = A_all(n,:);                   
        B_all_sink2(size(sink2,2),:) = B_all(n,:);                   
        
    end
end

if exist('A_all_sink1','var')
    avgA_sink1 = mean(A_all_sink1);                 avgB_sink1 = mean(B_all_sink1);
    sdevA_sink1 = std(A_all_sink1,0);               sdevB_sink1 = std(B_all_sink1,0);
end

if exist('A_all_sink2','var')
    avgA_sink2 = mean(A_all_sink2);                 avgB_sink2 = mean(B_all_sink2);
    sdevA_sink2 = std(A_all_sink2,0);               sdevB_sink2 = std(B_all_sink2,0);
end
%% Find unclassified trajectories
s = size(sink1,2) + size(sink2,2);

check = [];

if s ~= reps
    temp =  [sink1 sink2];
    ind = BubbleSort(temp);
    Dind = diff(ind);
    miss = find(Dind>1);
    for q=1:size(miss,2)
            check(q) = ind(miss(q))+1;
    end
end

check

%% Manually check why classification did not work.
% If variable 'check' is not empty, manually check indicated trajectories.
% Adjust lists 'sink1' and 'sink2' as needed.

s1 = size(sink1,2);                         s2 = size(sink2,2);
percON = s1 / (s1 + s2) * 100
percOFF = s2 / (s1 + s2) * 100

%% Graph ABK-stochastic and deterministic (selected) time trajectories
fig0 = figure('Name','Mutual Inhibition Time course','NumberTitle','off');
set(fig0,'Position',[501 1 500 406]);    
time = 0:dt:maxTime;
plot(time,A_all(2,:),'b');                                          hold on;
plot(time,B_all(2,:),'r');
plot(t_sol,y_sol(:,1),'--c');                  
plot(t_sol,y_sol(:,2),'--m');                               hold off;
axis([0 t_sol(end) 0 125]);                              % axis tight;
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('t (sec)');         ylabel('N(t)');         
leg0 = legend('ABK N_A(t)','ABK N_R(t)','DE N_A(t)','DE N_R(t)','Location','Best');
set(leg0,'Location','NorthWest');
set(leg0,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95]);

%% Plot State Space: nullclines, direction field, trajectories, separatrix, fixed points
fig1 = figure('Name','State Space','NumberTitle','off');                hold on;
set(fig1,'Position',[1 1 500 406]); 

% p_ncA = plot(A_nc,0:agents,'','Color',[1 0.3 0.2],...
%     'LineWidth',2,'DisplayName','N_A nullcline');
% p_ncR = plot(0:agents,B_nc,'','Color',[0.4 0.4 1],...
%     'LineWidth',2,'DisplayName','N_R nullcline');      % B = R

% streakarrow(Am,Bm,Ua,Ub,0.7,1);             % Direction field
% colorbar vert;             h=colorbar;    set(h,'ylim',[min(Vmag(:)) max(Vmag(:))]);
% quiver(Am,Bm,U,V,'g');                    % Direction field

% Plot separatrix
p_s = plot(0:agents,0:agents,'--','Color',[1 0.75 0],...
    'LineWidth',1.5,'DisplayName','Separatrix');

% Plot trajectories
p_de  = plot(y_sol(:,1),y_sol(:,2),'+k','DisplayName','DE  Trajectory','MarkerSize',5);  

if exist('A_all_sink1','var')
    p_abk1 = plot(avgA_sink1,avgB_sink1,'m','DisplayName','ABK Trajectory --> Sink 1');  
    p_abk1s0 = plot(avgA_sink1-sdevA_sink1,avgB_sink1-sdevB_sink1,':','Color',[0.8 0.8 0.8]);  
    p_abk1s1 = plot(avgA_sink1+sdevA_sink1,avgB_sink1+sdevB_sink1,':','Color',[0.8 0.8 0.8]);  
end

if exist('A_all_sink2','var')
    p_abk2 = plot(avgA_sink2,avgB_sink2,'c','DisplayName','ABK Trajectory --> Sink 2');  
    p_abk2s0 = plot(avgA_sink2-sdevA_sink2,avgB_sink2-sdevB_sink2,':','Color',[0.8 0.8 0.8]);  
    p_abk2s1 = plot(avgA_sink2+sdevA_sink2,avgB_sink2+sdevB_sink2,':','Color',[0.8 0.8 0.8]);  
end

if size(A_ss_int,1)==3
    for z=1:3
        if mod(z,2)==0              % Unstable fixed point
            p(z) = plot(A_ss_int(z),B_ss_int(z),'oc','MarkerSize',8,...
                'MarkerFaceColor','r','DisplayName','Unstable');
        else                        % Two Stable fixed points
            p(z) = plot(A_ss_int(z),B_ss_int(z),'oc','MarkerSize',8,...
                'MarkerFaceColor','b','DisplayName','Stable');
        end
    end
else
    plot(A_ss_int,B_ss_int,'oc','MarkerSize',8,...
        'MarkerFaceColor','b','DisplayName','Stable');          % Sole stable fixed point
end

xlabel('N_A');                    ylabel('N_R');                    
axis([0 110 0 110]);                                             hold off;  
set(gca,'XMinorTick','on','YMinorTick','on','Box','off','DataAspectRatio',[1 1 1]);

% leg1 = legend([p_ncA p_ncR p_s p(1) p(2)]);         % Without DE, ABK Trajectories

if exist('p_abk1','var') && exist('p_abk2','var')
    leg1 = legend([p_de p_abk1 p_abk2]);              % DE, ABK Trajectories
else, leg1 = legend();
end

set(leg1,'Location','NorthEast');
% set(leg1,'Position',[0.7 0.71 0.29 0.262]);
set(leg1,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95]);

%% Finish
close(pb);
clear n r w i j tempA tempB Av Bv P* p* leg* fig* T* temp z;
toc

%% Notes

% For Parameters:
% k_a = 1 = k_b                      % MICROSCOPIC basal rate of A and B synthesis (0th order)
% k_da = 0.01 = k_db                 % degradation rate constant (1st order) for A and B, resp.
% K_a = 50 = K_b                     % Populations producing half-maximal regulatory effect
% alpha_a = 0 = alpha_b              % Effect of A (B) on synthesis of B (A); (=0, Complete Repression)
% n_a = 3 = n_b                      % Hill Coefficient: measure of cooperativity

% sink 1: A=11.47 , B=98.9               sink 2: A=98.8 , B=11.47

%   Ai       Bi        -> sink1    -> sink2      # sim runs
%  ----     ----       --------    ---------     ----------
%   10       10         49.2 %      50.8 %          500
%   10       11         49.8 %      50.2 %          1000
%   10       90         99.975%     0.025%          4000   *
%   10       95         100  %      0    %          1000
%   10       100        100  %      0    %          2000
%   20       80         100  %      0    %          1000
%   30       75         99.8 %      0    %          500    (incl. unclassified run)
%   30       100        100  %      0    %          500
%   40       100        100  %      0    %          500
%   50       100        99.6 %      0.4  %          500
%   50       50         50.6 %      49.4 %          500
%   75       50         5.6  %      94.4 %          500

% * Starting from an initial condition close to sink 1 (Ai=10, Bi=90), the system goes to
% sink1 (A=11.47,B=98.9) 99.8% of the time, and to sink2 (A=98.8,B=11.47) 0.2% (n=500).
