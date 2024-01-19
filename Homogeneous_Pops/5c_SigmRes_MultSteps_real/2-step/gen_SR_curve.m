% This is the same as the Sigmoidal Response 2-step motif,
% except the steps occur serially, as would be expected in 
% a more realistic scenario.

% R <--> Rp <--> Rpp             Rates: k_f, k_r
%     ^       ^
%     |       |
%     S       S                  S: Signal (reactant in forward reactions)

% Assume number of S molecules does not change.
% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear;             tic;            clc; 
pb = waitbar(0,'0');

%% Declare variables
global k_f k_r S;

maxTime = 50;                         % Maximum Simulation time (sec)
dt = 1/100;                           % Constant (fixed) time step increment (sec)
t_steps = maxTime / dt;

agents = 100;                 	disp(['Agents = ' num2str(agents)]);
k_f = 0.01;                            % basal forward rate (2nd order)
k_r = 0.1;                             % basal reverse rate (1st order)

% ** Number of S molecules to do simulation for **
S_array = 0:40;
s = size(S_array,2);

% Preallocate memory for steady state values and rates
ABK_Rpp = zeros(1,s);

% Initialize - Preallocate memory for variables; ** Initial conditions **
Tr = zeros(s,t_steps);              % Total sum of R agents in each time step	
Trp = zeros(s,t_steps);             % Total sum of Rp agents in each time step
Trpp = zeros(s,t_steps);            % Total sum of Rpp agents in each time step

% ******** Initial conditions - Number of R, Rp, Rpp Agents ********
Tr(:,1) = floor(agents/4);       		Trp(:,1) = agents - Tr(1);
Trpp(1) = 0;
% ******************************************************************

%% Perform experiments for different S values
for n=1:s

    progress = n / s;
    waitbar(progress,pb,sprintf('%.0f%% %',progress*100));
    
    S = S_array(n);             % Number of S agents for this experiment

    tempR = zeros(1,agents);                          
    % Put "1" where agents are "alive", then randomize the array
    for c=1:Tr(1),                      tempR(c)=1;             end
    tempR = RandArray(tempR);           % Randomize array
    Rv = [tempR ; tempR];               % Initialize vector storing state of R agents     
    Rpv = ~Rv;                          % Initialize vector storing state of Rp agents
    Rppv = zeros(2,agents);             % Initialize vector storing state of Rpp agents
    % Notes on Rv, Rpv, Rppv:
    % - Markov process, so only previous and current time steps needed --> 2 rows:  
    % - R + Rp + Rpp = agents
    clear c tempR;

    % ***** ABK Simulation *****
    
    % ** Time-independent transition probabilities **
    P_f = 1 - exp(-k_f * S * dt);    % P_ber for forward reaction: 2nd order (but S is constant)            
%     P_f = k_f * S * dt;              % P_dif for forward reaction: 2nd order
    P_r = 1 - exp(-k_r * dt);        % P_ber for reverse reaction: 1st order            
%     P_r = k_r * dt;                  % P_dif for reverse reactio: 1st order

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
                if rand < P_r                  % **check probability condition**
                    Rppv(2,i) = 0;                 % Rpp agent "dies"
                    Rpv(2,i) = 1;                 % Rpp is converted to Rp
                end
            end
        end
        
        Tr(n,t+1) = sum(Rv(2,:));           Trp(n,t+1) = sum(Rpv(2,:));     
        Trpp(n,t+1) = sum(Rppv(2,:));
        
        Rv(1,:) = Rv(2,:);          Rpv(1,:) = Rpv(2,:);            Rppv(1,:) = Rppv(2,:); 
        
        t = t + 1;   
    end
    % Average Rp population over the last 25 sec of the simulation
    % Steady-state has been reached for this time period
    ABK_Rpp(n) = mean(Trpp(n,end-25*dt^-1:end));
    
end             % end "for n" loop

finaltime = (t-1) * dt;

%% Plot signal-response curve
fig1 = figure('Name','Signal-Response Curve','NumberTitle','off');
set(fig1,'Position',[1 1 500 406]);

% plot(S_array,Rp_ss,S_array,ABK_Rpp,'.r');
plot(S_array,ABK_Rpp,'.r');

% xlabel('S');    ylabel('R_{pp}');    
xlabel('N_S');    ylabel('N_R^*');              % For Report figure
axis([0 S_array(end) 0 agents]);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
% leg1 = legend('DE R_{p}^*','ABK R_{p}^*','Location','SouthEast');
% leg1 = legend('DE','ABK','Location','SouthEast');
% set(leg1,'Location','NorthWest');
% set(leg1,'FontName','Times New Roman','FontSize',9,...
%     'EdgeColor',[0.95 0.95 0.95]);

%% Find Coefficient of Determination, R^2
% % This works for any curve fit, linear or nonlinear.
% SST = sum((ABK_Rpp' - mean(ABK_Rpp)).^2);     % Total sum of squares for simulation data (for A)
% SSR = sum((ABK_Rpp' - Rp_ss').^2);           % sum of square residuals (sim data vs DE predictions)
% 
% Rsq = 1 - SSR./SST                          % Definition of R^2
% R = sqrt(Rsq);                              % Correlation Coefficient R
%% Determine Correlation between SR curve and theoretical expectation
% CC = corrcoef(Rp_ss,ABK_Rp);        % CC = CorrCoeff
% disp(['Corr Coeff R = ' num2str(CC(1,2))]);

%% Finish
close(pb);
clear pb progress fig* leg* p r w i j k n temp Rv Rpv;
toc
