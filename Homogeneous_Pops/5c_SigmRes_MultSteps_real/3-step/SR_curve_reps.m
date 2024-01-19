% This is the same as the Sigmoidal Response 2-step motif,
% except the steps occur serially, and one extra step has
% been added. This progression would be expected in 
% a more realistic scenario.

% R <--> Rp <--> Rpp <--> Rppp             Rates: k_f, k_r
%     ^       ^        ^
%     |       |        |
%     S       S        S                   S: Signal (reactant in forward reactions)

% Assume number of S molecules does not change.
% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear;             tic;            clc; 
rng(0);
pb = waitbar(0,'0');

%% Declare variables
global k_f k_r S;

maxTime = 100;                         % Maximum Simulation time (sec)
dt = 1/100;                           % Constant (fixed) time step increment (sec)
t_steps = maxTime / dt;

reps = 30;

agents = 100;                 	disp(['Agents = ' num2str(agents)]);
k_f = 0.01;                            % basal forward rate (2nd order)
k_r = 0.1;                             % basal reverse rate (1st order)

% ** Number of S molecules to do simulation for **
S_array = 0:40;
s = size(S_array,2);

% Preallocate memory for steady state values
ABK_R    = zeros(reps,s);           ABK_Rp   = zeros(reps,s);
ABK_Rpp  = zeros(reps,s);           ABK_Rppp = zeros(reps,s);

R_ss = zeros(1,s);                  Rp_ss = zeros(1,s);          
Rpp_ss = zeros(1,s);                Rppp_ss = zeros(1,s);
%% Perform experiments for different S values
for n=1:s

    progress = n / s;
    waitbar(progress,pb,sprintf('%.0f%% %',progress*100));
    
    % Initialize - Preallocate memory for variables; ** Initial conditions **
    Tr = zeros(reps,t_steps);              % Total sum of R agents in each time step	
    Trp = zeros(reps,t_steps);             % Total sum of Rp agents in each time step
    Trpp = zeros(reps,t_steps);            % Total sum of Rpp agents in each time step
    Trppp = zeros(reps,t_steps);           % Total sum of Rppp agents in each time step
    
    % ******** Initial conditions - Number of R, Rp, Rpp Agents ********
    Tr(:,1) = floor(agents/4);       		Trp(:,1) = agents - Tr(1);	
    Trpp(1) = 0;                            Trppp(1) = 0;
    % ******************************************************************

    S = S_array(n);             % Number of S agents for this experiment

    [t_sol, y_sol] = ode45(@sigm_3step_real_dif,0:0.2:maxTime,[Tr(1); Trp(1); Trpp(1); Trppp(1)]);
    R_ss(n)   = y_sol(end,1);               Rp_ss(n)  = y_sol(end,2); 
    Rpp_ss(n) = y_sol(end,3);               Rppp_ss(n) = y_sol(end,4);
    
    % ** Time-independent transition probabilities **
    P_f = 1 - exp(-k_f * S * dt);    % P_ber for forward reaction: 2nd order (but S is constant)            
    %     P_f = k_f * S * dt;              % P_dif for forward reaction: 2nd order
    P_r = 1 - exp(-k_r * dt);        % P_ber for reverse reaction: 1st order            
    %     P_r = k_r * dt;                  % P_dif for reverse reactio: 1st order

    for w=1:reps
        
        tempR = zeros(1,agents);                          
        % Put "1" where agents are "alive", then randomize the array
        for c=1:Tr(1),                      tempR(c)=1;             end
        tempR = RandArray(tempR);           % Randomize array
        Rv = [tempR ; tempR];               % Initialize vector storing state of R agents     
        Rpv = ~Rv;                          % Initialize vector storing state of Rp agents
        Rppv = zeros(2,agents);             % Initialize vector storing state of Rpp agents
        Rpppv = zeros(2,agents);            % Initialize vector storing state of Rppp agents
        
        % Notes on Rv, Rpv, Rppv, Rpppv:
        % - Markov process, so only previous and current time steps needed --> 2 rows:  
        % - R + Rp + Rpp + Rppp = agents
        clear c tempR;

        % ***** ABK Simulation *****
        t = 1;                              % Time counter variable

        while t <= t_steps

            for i = 1:agents
                if Rv(1,i) == 1                     % if R agent is alive
                    if rand < P_f                  % **check probability condition**
                        Rv(2,i) = 0;                 % R agent "dies"
                        Rpv(2,i) = 1;                 % R is converted to Rp
                    end

                elseif Rpv(1,i) == 1					% if Rp agent is "alive"
                    r = rand;
                    if r < P_r                  % **check probability condition**
                        Rpv(2,i) = 0;                 % Rp agent "dies"
                        Rv(2,i) = 1;                 % Rp is converted to R
                    elseif r >= P_r && r < P_r + P_f
                        Rpv(2,i) = 0;                 % Rp agent "dies"
                        Rppv(2,i) = 1;                 % Rp is converted to Rpp
                    end

                elseif Rppv(1,i) == 1					% if Rpp agent is "alive"
                    r = rand;
                    if r < P_r                  % **check probability condition**
                        Rppv(2,i) = 0;                 % Rpp agent "dies"
                        Rpv(2,i) = 1;                 % Rpp is converted to Rp
                    elseif r >= P_r && r < P_r + P_f
                        Rppv(2,i) = 0;                 % Rpp agent "dies"
                        Rpppv(2,i) = 1;                 % Rpp is converted to Rppp
                    end

                elseif Rpppv(1,i) == 1					% if Rppp agent is "alive"
                    if rand < P_r                  % **check probability condition**
                        Rpppv(2,i) = 0;                 % Rppp agent "dies"
                        Rppv(2,i) = 1;                 % Rppp is converted to Rpp
                    end
                end
            end            
            
            Tr(w,t+1) = sum(Rv(2,:));           Trp(w,t+1) = sum(Rpv(2,:));
            Trpp(w,t+1) = sum(Rppv(2,:));       Trppp(w,t+1) = sum(Rpppv(2,:));

            Rv(1,:) = Rv(2,:);                  Rpv(1,:) = Rpv(2,:);             
            Rppv(1,:) = Rppv(2,:);              Rpppv(1,:) = Rpppv(2,:);    
            
            t = t + 1;   
        end
        
        % Average populations over the last 50 sec of the simulation
        % Steady-state has been reached for this time period
        ABK_R(w,n)   = mean(Tr(w,end-50*dt^-1:end));
        ABK_Rp(w,n)  = mean(Trp(w,end-50*dt^-1:end));
        ABK_Rpp(w,n) = mean(Trpp(w,end-50*dt^-1:end));
        ABK_Rppp(w,n)= mean(Trppp(w,end-50*dt^-1:end));
        
    end             % end "for w" loop
end             % end "for n" loop

finaltime = (t-1) * dt;

avg_R   = mean(ABK_R);                        std_R = std(ABK_R);
avg_Rp  = mean(ABK_Rp);                       std_Rp = std(ABK_Rp);
avg_Rpp = mean(ABK_Rpp);                      std_Rpp = std(ABK_Rpp);
avg_Rppp= mean(ABK_Rppp);                     std_Rppp = std(ABK_Rppp);

%% Plot signal-response curve
fig1 = figure('Name','Signal-Response Curve','NumberTitle','off');
set(fig1,'Position',[1 1 500 406]);                             hold on;
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.

p_de = plot(S_array,Rppp_ss,'DisplayName','DE Triple-Step');
p_abk = plot(S_array,avg_Rppp,'.r','DisplayName','ABK Triple-Step Rppp');
plot(S_array,avg_Rppp+std_Rppp,'--','Color',[0.8 0.8 0.8]);
plot(S_array,avg_Rppp-std_Rppp,'--','Color',[0.8 0.8 0.8]);
% shadedErrorBar(S_array,avg_Rppp,std_Rppp,'r',1);

% plot(S_array,avg_R,'.b','DisplayName','ABK Dual-Step R');
% plot(S_array,avg_Rp,'xc','DisplayName','ABK Dual-Step Rp');
% plot(S_array,avg_Rpp,'xk','DisplayName','ABK Dual-Step Rpp');

% xlabel('S');    ylabel('R_ppp');    
xlabel('N_S');    ylabel('N_R^*');              % For Report figure Rppp = R

% Also plot hyperbolic curve if process wasn't multi-step, for comparison.
hyp = agents .* S_array ./ ( k_r/k_f + S_array );
p_hyp = plot(S_array,hyp,'g','DisplayName','DE Single-Step');

axis([0 S_array(end) 0 agents]);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
% leg1 = legend('DE R_{p}^*','ABK R_{p}^*','Location','SouthEast');
leg1 = legend([p_de , p_abk],'Location','SouthEast');
set(leg1,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95]);

%% Find Coefficient of Determination, R^2
% This works for any curve fit, linear or nonlinear.
SST = sum((avg_Rppp' - mean(avg_Rppp)).^2);     % Total sum of squares for simulation data (for A)
SSR = sum((avg_Rppp' - Rppp_ss').^2);           % sum of square residuals (sim data vs DE predictions)

Rsq = 1 - SSR./SST                          % Definition of R^2
% R = sqrt(Rsq);                              % Correlation Coefficient R
%% Determine Correlation between SR curve and theoretical expectation
% CC = corrcoef(Rp_ss,ABK_Rp);        % CC = CorrCoeff
% disp(['Corr Coeff R = ' num2str(CC(1,2))]);

%% Finish
close(pb);
clear pb progress fig* leg* p r w i j k n s S temp Rv Rpv Rppv p_* *_sol;
toc
