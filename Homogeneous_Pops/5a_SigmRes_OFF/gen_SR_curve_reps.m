% Sigmoidal Response system:                  [Note for Report: A <--> R]
% R <--> Rp             Rates: k_f, k_r
%    ^
%    |
%    S                 S: Signal (reactant in reverse reaction)
% Assume number of S molecules does NOT change.
% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;             tic;            clc; 
pb = waitbar(0,'0');

% Declare variables
global k_f k_r Km_f Km_r S;

reps = 10;

maxTime = 1000;                         % Maximum Simulation time (sec)
dt = 1/50;                              % Constant (fixed) time step increment (sec)
t_steps = maxTime / dt;

agents = 100;          	disp(['Agents = ' num2str(agents)]);
k_f = 1;                                % basal forward rate constant
k_r = 0.05;                             % basal reverse rate constant
Km_f = 5;                               % MICROSCOPIC MM constant for forward rx
Km_r = 30;                              % MICROSCOPIC MM constant for reverse rx

S_array = 1:60;                         % Number of S molecules is do simulation for
s = size(S_array,2);
Trp = zeros(s,t_steps,reps);            % Total sum of Rp agents in each time step

%% *** Set up differential equations of system using symbolic variables ***
for n=1:s
    S = S_array(n);
    syms R Rp positive;
    dR_sym =  - k_f * R / (Km_f + R)  +  k_r * Rp * S / (Km_r + Rp);
    dRp_sym = + k_f * R / (Km_f + R)  -  k_r * Rp * S / (Km_r + Rp);

    % Symbolic calculations - Calculate theoretical steady state values for R, Rp
    dR_sym_sub = subs(dR_sym,Rp,agents-R);       % Restriction: R+Rp=agents
    R_ss_sym = solve(dR_sym_sub == 0,R,'MaxDegree',4,'Real',true);
    R_ss = double(R_ss_sym);

    % Only want positive steady state values
    if size(R_ss,1)>1
        R_ss_all = R_ss;
        temp = find(R_ss <= agents & R_ss >= 0, 1);
        R_ss = R_ss_all(temp);
    end
    Rp_ss(n) = agents - R_ss;
end
        
%% Simulation
for a=1:reps
    
    progress = a / reps;
    waitbar(progress,pb,sprintf('%.0f%% %',progress*100));

    Tr = zeros(s,t_steps);              % Total sum of R agents in each time step	

    for n=1:s

        S = S_array(n);                         fprintf(1,'.');

    %  Initialize - Preallocate memory for variables; ** Initial conditions **
        t = 1;                              % Time counter variable
        P_f = zeros(1,t_steps);             % For storing probability value at each time step
        P_r = zeros(1,t_steps);             % For storing probability value at each time step

        % ******** Initial conditions - Number of R, Rp Agents ********
        Tr(:,1) = floor(agents/2);     		Trp(:,1,:) = agents - Tr(1);					
        % ************************************************************

        tempR = zeros(1,agents);                          
        % Put "1" where agents are "alive", then randomize the array
        for c=1:Tr(1,1),                    tempR(c)=1;             end
        tempR = RandArray(tempR);           % Randomize array
        Rv = [tempR ; tempR];               % Initialize vector storing state of R agents     
        Rpv = ~Rv;                          % Initialize vector storing state of Rp agents 
        % Notes on Rv, Rpv:
        % - Markov process, so only previous and current time steps needed --> 2 rows:          
        % - Rv and Rpv are complementary (R+Rp=agents)
        clear c tempR;

    %  ABK Simulation

        while t*dt <= maxTime               % Additional condition: abs(Tr(t)-R_ss) > 0.5

            P_f(t) = k_f / (Km_f + Tr(n,t)) * dt;           % P_dif forward rx
            P_r(t) = k_r * S / (Km_r + Trp(n,t,a)) * dt;    % P_dif reverse rx

            for i = 1:agents
                if Rv(1,i) == 1                             % if R agent is alive
                    if rand < P_f(t)                        % **check probability condition**
                        Rv(2,i) = 0;                        % R agent "dies"
                        Rpv(2,i) = 1;                       % R is converted to Rp
                    end
                elseif Rpv(1,i) == 1                        % if Rp agent is "alive"
                    if rand < P_r(t)                        % **check probability condition**
                        Rpv(2,i) = 0;                       % Rp agent "dies"
                        Rv(2,i) = 1;                        % Rp is converted to R
                    end
                end
            end

            Tr(n,t+1) = sum(Rv(2,:));                Trp(n,t+1,a) = sum(Rpv(2,:));

            Rv(1,:) = Rv(2,:);                       Rpv(1,:) = Rpv(2,:); 
            t = t + 1;   
        end

    %   Average Rp population over the last 400 sec of the simulation
    %   Steady-state has been reached for this time period
        

    end             % end "for n=1:s" loop
end                 % end "for a=1:reps" loop

avgRp = mean(Trp,3);                sdevRp = std(Trp,0,3);
ABK_Rp = mean(avgRp(:,end-400*dt^-1:end),2)';
ABK_Rp_dev = std(avgRp(:,end-400*dt^-1:end),0,2)';
%% Rate curve figure
fig0 = figure('Name','Rate Curve','NumberTitle','off');            
set(fig0,'Position',[550 1 500 406]);                    hold on;
Rp_array = 0:agents;
dRpdt_synt = k_f .* (agents - Rp_array) ./ (Km_f + (agents - Rp_array));

for n=1:s
    S = S_array(n); 
    
    dRpdt_degr(n,:) = k_r * S .* Rp_array ./ (Km_r + Rp_array); 

    [Rp_ss_int(n), ssr(n)] = intersections(Rp_array,dRpdt_degr(n,:),Rp_array,dRpdt_synt);
    % Rp_int: Steady-state value of Rp derived from intersections in Rate curve
    % ssr: Steady-state rate of forward and back reactions
    if mod(n,5) == 0
        p(n) = plot(Rp_array,dRpdt_degr(n,:),'-m','DisplayName','Degradation');
        plot(Rp_ss_int(n),ssr(n),'ok');
    end
end

p_s = plot(Rp_array,dRpdt_synt,'-','Color',[0 0.5 0],'DisplayName','Synthesis');
hold off;

xlabel('N_R');        ylabel('dN_R/dt');            % For report 
% xlabel('R_p');        ylabel('dR_p/dt');    
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
leg0 = legend([p(5) p_s],'Location','NorthWest');
set(leg0,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95]);

% Create textbox
annotation(fig0,'textbox',...
    [0.8 0.132 0.111332007952286 0.0519953917050691],...
    'String',{'N_S = 5'},...
    'FontWeight','demi','FontName','Times New Roman',...
    'FitBoxToText','off','LineStyle','none','Color',[0 0.75 0]);

% Create textbox
annotation(fig0,'textbox',...
    [0.8 0.48 0.121153081510935 0.0519953917050691],...
    'String',{'N_S = 30'},...
    'FontWeight','demi','FontName','Times New Roman',...
    'FitBoxToText','off','LineStyle','none','Color',[0 0.75 0]);

% Create textbox
annotation(fig0,'textbox',...
    [0.8 0.86 0.121153081510935 0.0519953917050691],'String',...
    {'N_S = 60'},...
    'FontWeight','demi','FontName','Times New Roman',...
    'FitBoxToText','off','LineStyle','none','Color',[0 0.75 0]);

%% Plot signal-response curve
fig1 = figure('Name','Signal-Response Curve','NumberTitle','off');
set(fig1,'Position',[1 1 500 406]);                    hold on;
p1 = plot(S_array,Rp_ss,'--b','DisplayName','DE');
plot(S_array,ABK_Rp,'.r');
% errorbar(S_array,ABK_Rp,ABK_Rp_dev,'.r');         % Using Matlab built-in function
% p2 = shadedErrorBar(S_array,ABK_Rp,ABK_Rp_dev,'r');    % Using File Exchange script

xlabel('N_S');                              % ylabel('N_{Rp}^*');  
ylabel('N_R^*');                            % For Report    
axis([0 S_array(end) 0 agents]);                        hold off;
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
leg1 = legend('DE','ABK');
set(leg1,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');

%% Find Coefficient of Determination, R^2
% This works for any curve fit, linear or nonlinear.
SST = sum((ABK_Rp' - mean(ABK_Rp)).^2);     % Total sum of squares for simulation data (for Rp)
SSR = sum((ABK_Rp' - Rp_ss').^2);           % sum of square residuals (sim data vs DE predictions)

Rsq = 1 - SSR./SST                         % Definition of R^2
% R = sqrt(Rsq);                             % Correlation Coefficient R

%% SR Deviation plot

deviation = ABK_Rp - Rp_ss;
figure; plot(deviation,'.r');
xlabel('N_S');                              
ylabel('ABK N_R^* - DE N_R^*');                            % For Report    

large_dev = find(abs(deviation)>5);

%% Finish
close(pb);
clear progress a r w i j k n t S s temp R_ss R_ss_all Rv Rpv P_f P_r p* fig* leg*;
toc
