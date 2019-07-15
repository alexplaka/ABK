% This is the same as the Hyperbolic Response system,
% except the forward rx depends on S^2, as would be 
% the case if two molecules of S were needed to 
% accomplish the conversion R --> Rp.

% R <--> Rp             Rates: k_f, k_r
%     ^
%     |
%     2S                S: Signal (reactant in forward reaction)
% Assume number of S molecules does NOT change.
% Agent-based interpretation of S population size.
% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;             tic;            clc; 
pb = waitbar(0,'0');

%% Declare variables
maxTime = 50;                         % Maximum Simulation time (sec)
dt = 1/100;                           % Constant (fixed) time step increment (sec)
t_steps = maxTime / dt;

agents = 100;                 	disp(['Agents = ' num2str(agents)]);
k_f = 0.0001;                            % basal forward rate (2nd order)
k_r = 0.1;                             % basal reverse rate (1st order)

% ** Number of S molecules to do simulation for **
S_array = 0:80;
s = size(S_array,2);

% Preallocate memory for steady state values and rates
ABK_Rp = zeros(1,s);

% Initialize - Preallocate memory for variables; ** Initial conditions **
Tr = zeros(s,t_steps);              % Total sum of R agents in each time step	
Trp = zeros(s,t_steps);             % Total sum of Rp agents in each time step

% ******** Initial conditions - Number of R, Rp Agents ********
Tr(:,1) = floor(agents/4);       		Trp(:,1) = agents - Tr(1);					
% ************************************************************

%% Perform experiments for different S values
for n=1:s

    progress = n / s;
    waitbar(progress,pb,sprintf('%.0f%% %',progress*100));
    
    S = S_array(n);             % Number of S agents for this experiment

    tempR = zeros(1,agents);                          
    % Put "1" where agenst are "alive", then randomize the array
    for c=1:Tr(1,1),                     tempR(c)=1;             end
    tempR = RandArray(tempR);           % Randomize array
    Rv = [tempR ; tempR];               % Initialize vector storing state of R agents     
    Rpv = ~Rv;                          % Initialize vector storing state of Rp agents 
    % Notes on Rv, Rpv:
    % - Markov process, so only previous and current time steps needed --> 2 rows:          
    % - Rv and Rpv are complementary (R+Rp=agents)
    clear c tempR;

    % ***** ABK Simulation *****
    
    % ** Time-independent transition probabilities **
    P_f = 1 - exp(-k_f * S*(S-1) * dt);    % P_ber for forward reaction: 1st order            
    %     P_f = k_f * S*(S-1) * dt;              % P_dif for forward reaction: 1st order
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
                if rand < P_r                   % **check probability condition**
                    Rpv(2,i) = 0;                 % Rp agent "dies"
                    Rv(2,i) = 1;                 % Rp is converted to R
                end
            end
        end
        
        Tr(n,t+1) = sum(Rv(2,:));                
        Trp(n,t+1) = sum(Rpv(2,:));
        
        Rv(1,:) = Rv(2,:);                       Rpv(1,:) = Rpv(2,:); 
        t = t + 1;   
    end
    % Average Rp population over the last 25 sec of the simulation
    % Steady-state has been reached for this time period
    ABK_Rp(n) = mean(Trp(n,end-25*dt^-1:end));
    
end             % end "for n" loop

finaltime = (t-1) * dt;
%% Rate Curve figure
fig0 = figure('Name','Rate Curve','NumberTitle','off');            
set(fig0,'Position',[550 1 500 406]);                    hold on;

Rp_array = 0:agents;
dRpdt_synt = zeros(s,size(Rp_array,2));
Rp_ss = zeros(1,s);                 ssr = zeros(1,s);

for n=1:s
    S = S_array(n);             % Number of S agents for this experiment
    
    Rp_ss(n) = agents * S*(S-1) / ( k_r/k_f + S*(S-1) );
    ssr(n) = k_r * Rp_ss(n);
    % ssr: Steady-state rate of forward and back reactions
    
    dRpdt_synt(n,:) = k_f * S*(S-1) .* (agents - Rp_array);  

    if mod(S,5) == 0
        p(n) = plot(Rp_array,dRpdt_synt(n,:),'-b','DisplayName','Synthesis');
    end
end

% Degradation curve is the same for all S
dRpdt_degr = k_r .* Rp_array;
p_d = plot(Rp_array,dRpdt_degr,'r');  
set(p_d,'DisplayName','Degradation');

for n=1:s
    S = S_array(n);
    if mod(S,5) == 0        
        scatter(Rp_ss(n),ssr(n),'ok');    
    end
end

hold off;
% xlabel('R_p');        ylabel('dR_p/dt');  
xlabel('N_R');        ylabel('dN_R / dt');              % For Report figure
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
leg0 = legend([p(1) p_d]);
set(leg0,'Location','NorthEast');
set(leg0,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95]);

% Create textbox
annotation(fig0,'textbox',...
    [0.121272365805168 0.209677419354839 0.111332007952286 0.0519953917050691],...
    'String',{'N_S = 5'},...
    'FontWeight','demi','FontName','Times New Roman',...
    'FitBoxToText','off','LineStyle','none','Color',[0 0.5 0]);

% Create textbox
annotation(fig0,'textbox',...
    [0.121391650099403 0.473087557603686 0.121153081510935 0.0519953917050691],...
    'String',{'N_S = 20'},...
    'FontWeight','demi','FontName','Times New Roman',...
    'FitBoxToText','off','LineStyle','none','Color',[0 0.5 0]);

% Create textbox
annotation(fig0,'textbox',...
    [0.123499005964213 0.900092165898615 0.121153081510935 0.0519953917050691],'String',...
    {'N_S = 40'},...
    'FontWeight','demi','FontName','Times New Roman',...
    'FitBoxToText','off','LineStyle','none','Color',[0 0.5 0]);

%% Plot signal-response curve
fig1 = figure('Name','Signal-Response Curve','NumberTitle','off');
set(fig1,'Position',[1 1 500 406]);
plot(S_array,Rp_ss,S_array,ABK_Rp,'.r');
% xlabel('S');    ylabel('R_p');    
xlabel('N_S');    ylabel('N_R^*');              % For Report figure
axis([0 S_array(end) 0 agents]);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
% leg1 = legend('DE R_{p}^*','ABM R_{p}^*','Location','SouthEast');
leg1 = legend('DE','ABK','Location','SouthEast');
set(leg1,'Location','NorthWest');
set(leg1,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95]);

%% Find Coefficient of Determination, R^2
% This works for any curve fit, linear or nonlinear.
SST = sum((ABK_Rp' - mean(ABK_Rp)).^2);     % Total sum of squares for simulation data (for A)
SSR = sum((ABK_Rp' - Rp_ss').^2);           % sum of square residuals (sim data vs DE predictions)

Rsq = 1 - SSR./SST                          % Definition of R^2
% R = sqrt(Rsq);                              % Correlation Coefficient R
%% Determine Correlation between SR curve and theoretical expectation
% CC = corrcoef(Rp_ss,ABK_Rp);        % CC = CorrCoeff
% disp(['Corr Coeff R = ' num2str(CC(1,2))]);

%% Finish
close(pb);
clear pb progress fig* leg* p r w i j k n temp Rv Rpv;
toc
