% This is the same as the Hyperbolic Response system,
% except the forward rx depends on S^2, as would be 
% the case if two molecules of S were needed to 
% accomplish the conversion R --> Rp.

% R <--> Rp             Rates: k_f, k_r
%     ^
%     |
%     2S                S: Signal (reactant in forward reaction)
% Assume number of S molecules does NOT change.
% 
% Here, we consider compositional heterogeneity in k_f or k_r (or both).

% Here, we construct the SR curve for this motif.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;             tic;            clc; 
pb = waitbar(0,'0');

%% Declare variables
maxTime = 100;                         % Maximum Simulation time (sec)
dt = 1/100;                           % Constant (fixed) time step increment (sec)
t_steps = maxTime / dt;

reps = 100;

agents = 100;                 	disp(['Agents = ' num2str(agents)]);

k_f_mean = 0.0001;                     % mean basal forward rate (2nd order)
k_f_sdev = 0.000025;                   % sdev of basal forward rate (2nd order)
% Set up R population structure (Choose one of the below 3 options):
k_f = k_f_mean * ones(1,agents);                % Homogeneous Population
% k_f = k_f_mean + k_f_sdev * randn(1,agents);  % PHM for Gaussian k_f values
% k_f = [k_f_mean/4 * ones(1,agents/2) , k_f_mean*7/4 * ones(1,agents/2)]; % PHM for preset k_f values


k_r_mean = 0.1;                        % mean basal reverse rate (1st order)
k_r_sdev = 0.05;                       % sdev of basal reverse rate (1st order)
% Set up Rp population structure (Choose one of the below 3 options):
% k_r = k_r_mean * ones(1,agents);                % Homogeneous Population
% k_r = [k_r_mean/4 * ones(1,agents/2) , k_r_mean*7/4 * ones(1,agents/2)]; % PHM for preset k_r values
k_r = k_r_mean + k_r_sdev * randn(1,agents);  % PHM for Gaussian k_r values
while isempty(find(k_r<0, 1)) == false      % Make sure values in PHM are >0 
    k_r = k_r_mean + k_r_sdev * randn(1,agents);
end

het_details_kf = het_Stats(k_f);
het_details_kr = het_Stats(k_r);

% ** Number of S molecules to do simulation for **
S_array = 0:80;
s = size(S_array,2);

% Preallocate memory for steady state values
Rp_ss = zeros(1,s);
ABK_Rp = zeros(reps,s);
avg_k_f = zeros(s,t_steps+1);       % For storing instantaneous <k_f(t)> for all values of S
avg_k_r = zeros(s,t_steps+1);       % For storing instantaneous <k_r(t)> for all values of S

%% Perform experiments for different S values
for n=1:s

    progress = n / s;
    waitbar(progress,pb,sprintf('%.0f%% %',progress*100));
        
    % Initialize - Preallocate memory for variables; ** Initial conditions **
    Tr = zeros(reps,t_steps);              % Total sum of R agents in each time step	
    Trp = zeros(reps,t_steps);             % Total sum of Rp agents in each time step

    % ******** Initial conditions - Number of R, Rp Agents ********
    Tr(:,1) = floor(agents/4);       		Trp(:,1) = agents - Tr(1);					
    % ************************************************************

    avg_k_f_temp = zeros(reps,t_steps+1);  % For storing instantaneous <k_f(t)> in 'reps' loop
    avg_k_r_temp = zeros(reps,t_steps+1);  % For storing instantaneous <k_r(t)> in 'reps' loop
    
    S = S_array(n);             % Number of S agents for this experiment

    Rp_ss(n) = agents * S^2 / ( k_r_mean/k_f_mean + S^2 );   % Deterministic steady-state

    % ** Time-independent transition probabilities **
    P_f = 1 - exp(-k_f * S^2 * dt);    % P_ber for forward reaction: 1st order            
    %     P_f = k_f * S^2 * dt;              % P_dif for forward reaction: 1st order
    P_r = 1 - exp(-k_r * dt);        % P_ber for reverse reaction: 1st order (using k_r PHM)            
    %     P_r = k_r * dt;                  % P_dif for reverse reactio: 1st order (using k_r PHM) 

    for w=1:reps
        
        tempR = zeros(1,agents);                          
        % Put "1" where agenst are "alive", then randomize the array
        for c=1:Tr(1,1),                     tempR(c)=1;             end
        tempR = RandArray(tempR);           % Randomize array
        Rv = [tempR ; tempR];               % Initialize vector storing state of R agents     
        Rpv = ~Rv;                          % Initialize vector storing state of Rp agents 
        % Notes on Rv, Rpv:
        % - Markov process, so only previous and current time steps needed --> 2 rows:          
        % - Rv and Rpv are complementary (R+Rp=agents)
        
        tempR = Rv(2,:)==1;           % Logical array of initial existing R agents
        avg_k_f_temp(w,1) = tempR * k_f' / sum(tempR);
        
        tempRp = Rpv(2,:)==1;         % Logical array of initial existing Rp agents
        avg_k_r_temp(w,1) = tempRp * k_r' / sum(tempRp); 

        clear c tempR*;

        % ***** ABK Simulation *****
        t = 1;                              % Time counter variable

        while t <= t_steps

            for i = 1:agents
                if Rv(1,i) == 1                     % if R agent is alive
                    if rand < P_f(i)                  % **check probability condition**
                        Rv(2,i) = 0;                 % R agent "dies"
                        Rpv(2,i) = 1;                 % R is converted to Rp
                    end
                elseif Rpv(1,i) == 1					% if Rp agent is "alive"
                    if rand < P_r(i)                % **check probability condition**
                        Rpv(2,i) = 0;                 % Rp agent "dies"
                        Rv(2,i) = 1;                 % Rp is converted to R
                    end
                end
            end

            Tr(w,t+1) = sum(Rv(2,:));                
            Trp(w,t+1) = sum(Rpv(2,:));

            Rv(1,:) = Rv(2,:);                       Rpv(1,:) = Rpv(2,:); 
            
            tempR = Rv(2,:)==1;     % Logical array of existing R agents at end of time step
            avg_k_f_temp(w,t+1) = tempR * k_f' / sum(tempR); % Find <k_f> at end of time step
            
            tempRp = Rpv(2,:)==1;     % Logical array of existing Rp agents at end of time step
            avg_k_r_temp(w,t+1) = tempRp * k_r' / sum(tempRp); % Find <k_r> at end of time step

            t = t + 1;   
        end
        
        avg_k_f(n,:) = mean(avg_k_f_temp,'omitNaN');    % <k_f> for each sampled S value
        avg_k_r(n,:) = mean(avg_k_r_temp,'omitNaN');    % <k_r> for each sampled S value
        
        % Average Rp population over the last 50 sec of the simulation
        % Steady-state has been reached for this time period
        ABK_Rp(w,n) = mean(Trp(w,end-50*dt^-1:end));
        
    end             % end "for w" loop
end             % end "for n" loop

finaltime = (t-1) * dt;
time = 0:dt:finaltime;

avg_Rp = mean(ABK_Rp);                      sdev_Rp = std(ABK_Rp);

%% Plot signal-response curve
fig1 = figure('Name','Signal-Response Curve','NumberTitle','off');
set(fig1,'Position',[1 1 500 406]);                             hold on;
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.

p_de = plot(S_array,Rp_ss,'DisplayName','HOM');
p_abk = plot(S_array,avg_Rp,'.r','DisplayName','HET');
plot(S_array,avg_Rp+sdev_Rp,'--','Color',[0.8 0.8 0.8]);
plot(S_array,avg_Rp-sdev_Rp,'--','Color',[0.8 0.8 0.8]);
% shadedErrorBar(S_array,avg_Rp,std_Rp,'r',1);

% ylabel('R_p');    
xlabel('N_S');    ylabel('N_R^*');              % For Report figure Rp = R

% % Also plot hyperbolic curve if process wasn't multi-step, for comparison.
% hyp = agents .* S_array.^1 ./ ( k_r/k_f + S_array.^1 );
% p_hyp = plot(S_array,hyp,'g','DisplayName','DE Single-Step');

axis([0 S_array(end) 0 agents]);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
% leg1 = legend('DE R_{p}^*','ABK R_{p}^*','Location','SouthEast');
leg1 = legend([p_de , p_abk],'Location','SouthEast');
set(leg1,'FontName','Times New Roman','FontSize',10,...
    'EdgeColor',[0.95 0.95 0.95]);

%% Plot <k_r> for a range of S values
fig2 = figure('Name','SigmoidalSwitch: avg k_r','NumberTitle','off');
set(fig2,'Position',[1 1 500 406]);                             hold on;
fig2.PaperUnits = 'inches';
fig2.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.

plot(time(1:10:end),avg_k_r(2:5:end,1:10:end)*10^2);

axis([time(1) time(end) 2 11.2]);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('$t \, \textrm{(sec)}$','Interpreter','LateX','FontSize',11);         
ylabel('$< k_r > \times 10^2$','Interpreter','LateX','FontSize',11);    

% Create arrow
annotation(fig2,'arrow',[0.92 0.92],[0.15 0.75]);

% Create textboxes
annotation(fig2,'textbox',[0.89 0.105 0.11 0.05],'String','$N_S=0$',...
    'Interpreter','latex','FitBoxToText','off','EdgeColor','none');

annotation(fig2,'textbox',[0.89 0.75 0.11 0.05],'String','$N_S=80$',...
    'Interpreter','latex','FitBoxToText','off','EdgeColor','none');

%% Plot <k_f> for a range of S values
fig3 = figure('Name','SigmoidalSwitch: avg k_f','NumberTitle','off');
set(fig3,'Position',[1 1 500 406]);                             hold on;
fig3.PaperUnits = 'inches';
fig3.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.

plot(time(1:10:end),avg_k_f(2:5:end,1:10:end)*10^5);

axis([time(1) time(end) 2 11.2]);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('$t \, \textrm{(sec)}$','Interpreter','LateX','FontSize',11);         
ylabel('$< k_f > \times 10^5$','Interpreter','LateX','FontSize',11);    

% Create arrow
annotation(fig3,'arrow',[0.92 0.92],[0.83 0.39]);

% Create textboxes
annotation(fig3,'textbox',[0.89 0.34 0.11 0.05],'String','$N_S=80$',...
    'Interpreter','latex','FitBoxToText','off','EdgeColor','none');

annotation(fig3,'textbox',[0.89 0.83 0.11 0.05],'String','$N_S=0$',...
    'Interpreter','latex','FitBoxToText','off','EdgeColor','none');


%% Find Coefficient of Determination, R^2
% This works for any curve fit, linear or nonlinear.
Rsq = CoefDet(avg_Rp,Rp_ss);                   
% R = sqrt(Rsq);                              % Correlation Coefficient R
%% Determine Correlation between SR curve and theoretical expectation
% CC = corrcoef(Rp_ss,ABK_Rp);        % CC = CorrCoeff
% disp(['Corr Coeff R = ' num2str(CC(1,2))]);

%% Finish
close(pb);
clear pb progress fig* leg* p p_* r w i j k n s S t *temp* Rv Rpv Tr Trp;
toc
