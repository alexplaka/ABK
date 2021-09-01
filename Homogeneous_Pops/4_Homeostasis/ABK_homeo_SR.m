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

clear;      clc;     tic;                           
rng(1);

global maxRagents agents k_b k_d k_f k_r Km_f Km_r S;

totalTime = 1500;                   % Simulation time (sec)
dt = 1/50;                          % Constant (fixed) time step increment (sec)
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
S_array = 1:40;                    % Do simulation for this range of S values

% Preallocate memory
R_ss = zeros(1,size(S_array,2));            E_ss = zeros(1,size(S_array,2));
Rv = zeros(size(S_array,2),t_steps);        Ev = zeros(size(S_array,2),t_steps);
R_end = zeros(1,size(S_array,2));           

for n = 1:size(S_array,2)
    
    S = S_array(n);                    % Assume number of S molecules/agents is NOT changing
    disp(['S=' num2str(S)]);
    
%     [temp1, temp2] = RateCurve;
%     % For S=0, there's no steady-state value. R increases monotonically.
%     if isempty(temp1)==1  
%         temp1 = maxRagents;        temp2 = NaN;
%     end
%     R_ss(n) = temp1;                    E_ss(n) = temp2;

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

    t = 1;

    P_d = k_d * S * dt;          % Probability of R degradation wrt R

    while t < t_steps % && Tr(t) < maxRagents
        P_b(t) = k_b * Te(t) * dt;   % Probability of R synthesis: 0th order UPreg by E
        P_f(t) = k_f / (Km_f + Tep(t)) * dt;          % wrt each Ep molecule
        P_r(t) = k_r * Tr(t) / (Km_r + Te(t)) * dt;   % wrt each E molecule

        % Take care of 0th order processes first
    %   treat R synthesis as 0th order, UPregulated by E (b/c E is not consummed)    
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
    
    R_end(n) = mean(Tr(end-750/dt:end));    % Average population size during last 750 sec
    
    % Remove unnecessary terminal 0's from arrays
%     if t < t_steps
%         Tr = Tr(1:t);                   Tr(t+1:t_steps) = NaN;
%         Tep = Tep(1:t);                 Tep(t+1:t_steps) = NaN;
%         Te = Te(1:t);                   Te(t+1:t_steps) = NaN;    
%         P_f = P_f(1:t);        P_r = P_r(1:t);         P_b = P_b(1:t); 
%     end

    Rv(n,:) = Tr;                           Ev(n,:) = Te;

end                 % end 'for n'

finaltime = (t-1) * dt;
[t_sol, y_sol] = ode45(@homeo_dif,0:finaltime/500:finaltime,[Tr(1); Tep(1); Te(1)]);

%% Rate Curve figure
dRdt_degr = zeros(size(S_array,2),maxRagents);          % Preallocate memory
ssr = zeros(1,size(S_array,2));                         % Preallocate memory
p_d = zeros(1,size(S_array,2));

fig0 = figure('Name','Rate Curve','NumberTitle','off');            
set(fig0,'Position',[550 1 500 406]);                    hold on;

syms E Ep R positive;
dEp_sym = + k_r * E * R / (Km_r + E) - k_f * Ep / (Km_f + Ep);
temp = subs(dEp_sym,E,agents-Ep);
Ep_sym = solve(temp==0,Ep);                 % result has 2 entries, -/+
% Denominator contains R - k_f/k_r

R_array = 0:maxRagents;    
k = k_f / k_r;
if k < maxRagents 
    R_array(k_f/k_r + 1) = [];      % Remove R=k_f/k_r entry (avoid div by 0)
end

Ep_array = double(subs(Ep_sym(2),R,R_array));
E_array = agents - Ep_array;

dRdt_synt = k_b .* E_array; 

for n=1:size(S_array,2)
    S = S_array(n);             % Number of S agents for this experiment
    
    dRdt_degr(n,:) = k_d .* R_array * S;

    [R_ss(n), ssr(n)] = intersections(R_array,dRdt_degr(n,:),R_array,dRdt_synt);
    
    if mod(S,5) == 0
        p_d(n) = plot(R_array,dRdt_degr(n,:),'m','DisplayName','Degradation');  
    end
end

% Synthesis curve is the same for all S values
p_s = plot(R_array,dRdt_synt,'-','Color',[0 0.5 0],'DisplayName','Synthesis');

for n=1:size(S_array,2)
    S = S_array(n);
    if mod(S,5) == 0        
        scatter(R_ss(n),ssr(n),'ok');    
    end
end

hold off;
xlabel('N_R');        ylabel('dN_R / dt');              % For Report figure
axis([0 50 0 2.5]);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
leg0 = legend([p_s p_d(5)]);
set(leg0,'Location','NorthWest');
set(leg0,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95]);

% Create textbox
annotation(fig0,'textbox',...
    [0.90 0.25 0.111332007952286 0.0519953917050691],...
    'String',{'N_S = 5'},...
    'FontWeight','demi','FontName','Times New Roman',...
    'FitBoxToText','off','LineStyle','none','Color',[0 0.75 0]);

% Create textbox
annotation(fig0,'textbox',...
    [0.90 0.74 0.121153081510935 0.0519953917050691],...
    'String',{'N_S = 20'},...
    'FontWeight','demi','FontName','Times New Roman',...
    'FitBoxToText','off','LineStyle','none','Color',[0 0.75 0]);

% Create textbox
annotation(fig0,'textbox',...
    [0.56 0.92 0.121153081510935 0.0519953917050691],'String',...
    {'N_S = 40'},...
    'FontWeight','demi','FontName','Times New Roman',...
    'FitBoxToText','off','LineStyle','none','Color',[0 0.75 0]);

%% Plot (selected) time trajectories
% time = 0:dt:finaltime;
% scatter(time,Rv(21,:),3,'.b');                                               
% scatter(time,Ev(1,:),3,'og');
% scatter(t_sol,y_sol(:,1),3,'xr');                                       hold on;
% scatter(t_sol,y_sol(:,3),3,'.g'); 
% axis([0 finaltime 0 agents]);         xlabel('time');            ylabel('#');     hold off; 
% legend('R stoc','E stoc','R deter','E deter','Location','Best');

%% Plot State Space - Modify accordingly
% figure('Name','State Space','NumberTitle','off');               hold on;
% scatter(y_sol(:,1),y_sol(:,3),'+k');                 % plot deterministic trajectory
% 
% % Modify the following for plotting specific trajectories
% comet(Rv(29,:),Ev(29,:));                    % plot stochastic trajectory
% 
% plot(R_ss,E_ss,'.c','MarkerSize',6);
% xlabel('R');                    ylabel('E');                  hold off; 
% legend('Deterministic','ABK stochastic');
%% Plot SR Curve
fig1 = figure('Name','Signal-Response Curve','NumberTitle','off');
set(fig1,'Position',[1 1 500 406]);
plot(S_array,R_ss,S_array,R_end,'.r');
xlabel('N_S');    ylabel('N_R^*');              
axis([1 S_array(end) 0 70]);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
leg1 = legend('DE','ABK');
set(leg1,'Location','NorthEast');
set(leg1,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95]);

%% Find Coefficient of Determination, R^2
% This works for any curve fit, linear or nonlinear.
SST = sum((R_end' - mean(R_end)).^2);      % Total sum of squares for simulation data (for A)
SSR = sum((R_end' - R_ss').^2);            % sum of square residuals (sim data vs DE predictions)

Rsq = 1 - SSR./SST                         % Definition of R^2
% R = sqrt(Rsq);                             % Correlation Coefficient R

%% Finish - Notes
clear h i j k n p* t Tr Tep Te R Ep E fig* leg*;
clear temp1 temp2 temp tempR tempEp tempE;
toc
