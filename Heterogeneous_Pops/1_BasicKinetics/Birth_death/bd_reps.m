%  --> A -->  
%  k_b   k_d
% Simulating birth-death process for a ** heterogeneous ** population of A: 
% k_b: birth, 0th order process 
% k_d: death, 1st order process: normally-distributed. 

% WARNING: This script can generate more than 8GB of data depending on the 
% chosen parameters (e.g., reps). Make sure you have enough RAM and swap
% so you don't get an 'out of memory' error message.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear;            tic;
rng(0);
pb = waitbar(0,'0');
%% Declare variables and functions
% ******** Initial condition - Number of A agents ********
Ao = 0;                       % Initial population size of A

reps = 1000;                  % Repeat simulation this many times

maxTime = 500;                % Maximum Simulation time (sec)  ********************
dt = 1/100;                   % Constant (fixed) time step increment (sec)
t_steps = maxTime / dt;

k_b = 0.1;                    % Microscopic 0th order birth rate (units: 1/sec)
k_d_mean = 0.01;              % Average "Death" rate; 1st order (units: 1/sec)
k_d_std = 0.0025;             % SDev of "Death" rate

A_ss = k_b / k_d_mean;        % Steady-state value
disp(['A_ss = ' num2str(A_ss)]);

agents = floor(4 * A_ss);

% *** Preallocate memory ***
k_d_instavg_all = zeros(reps,t_steps);  % "instanteneous" mean k_d values
k_d_inststd_all = zeros(reps,t_steps);  % sdev of "instanteneous" k_d values

Sa = zeros(reps,t_steps);     
% **************************

% ***** Probability of birth condition can be specified here ******
P_b = k_b * dt;               % Probability of "birth" process (0th order)
% *********************************************************************************

%% ABK Simulation 

for n=1:reps
    
    progress = n/reps;
    waitbar(progress,pb,sprintf('%.0f%% %',progress*100));

    % Initialize - Preallocate memory for variables; ** Initial conditions **
    At = zeros(1,t_steps);              % Sum of A agents in each time step	

    At(1) = Ao;        
    
    k_d = zeros(1,agents);
    k_d_all = zeros(t_steps,agents);  % for collecting agent-specific k_d values
    % ************************************************************

    tempA = zeros(1,agents);            % tempB = zeros(1,agents);              
    % Put "1" where agents are "alive", then randomize the array
    for c=1:At(1),                      tempA(c)=1;             end
    tempA = RandArray(tempA);           % Randomize array
    Av = [tempA ; tempA];               % Initialize vector storing state of A agents     
    % - Markov process, so only previous and current time steps needed --> 2 rows:          
    clear c tempA;

    time = zeros(1,t_steps);
    
    for t=2:t_steps
    %     dt = 1 / (ko*Sa(t-1));            % Variable time step increment ; OLD
    %     dt = exprnd(1/(k*S(t-1)));   % Exponentially-distributed time step increment ; OLD

        if rand < P_b                           % "Birth", 0th order reaction
            temp = find(Av(1,:)==0);            % Randomly choose agent which becomes "alive" 
            ch = temp(ceil(rand * size(temp,2)));
            Av(2,ch) = 1;
            k_d(ch) = k_d_mean + k_d_std * randn;  
            while k_d(ch) < 0
                k_d(ch) = k_d_mean + k_d_std * randn;  % Retry
                disp(['Negative k_d! --> New k_d = ' num2str(k_d(ch))]);  
            end
        end

        tempA1 = find(Av(1,:)==1);
        for i = 1:size(tempA1,2)
%             P_d = k_d(tempA1(i)) * dt;               % Probability of "death" process (1st order), P_dif
            P_d = 1 - exp(-k_d(tempA1(i)) * dt);     % Probability of "death" process (1st order), P_ber

            if rand < P_d                       % "Death", 1st order reaction
                Av(2,tempA1(i)) = 0;            % A agent dies
                k_d(tempA1(i)) = 0;
            end
        end

        At(t) = sum(Av(2,:)); 
        k_d_all(t,:) = k_d;
        
        Av(1,:) = Av(2,:);                     
        time(t) = time(t-1) + dt;
        
    end
    
    Sa(n,:) = At; 
    
    for x=1:agents
        for y=1:size(time,2)
            if k_d_all(y,x) == 0
                k_d_all(y,x) = NaN;
            end
        end
    end
    
    k_d_instavg_all(n,:) = mean(k_d_all,2,'omitNaN');
    k_d_inststd_all(n,:) = std(k_d_all,0,2,'omitNaN');
    
end         % end "for n" loop

k_d_instavg  = mean(k_d_instavg_all,'omitNaN');
k_d_inststd = mean(k_d_inststd_all,'omitNaN');

%% See if extinction (N_A(t)=0) occurred in any of the trajectories
for w=1:size(Sa,1)
    if Sa(w,end)==0
        disp(['rep = ' num2str(w) ' ; Extinction!']);
    end
end

clear k_d k_d_all i n At Av ch t temp tempA1 w x y;
%% Calculate Average and SDev of time trajectory
avgA = mean(Sa);                     sdevA = std(Sa);               
cvA = sdevA ./ avgA;                                % Coefficient of variation

%% Compute deterministic trajectory (using analytic solution)
t_sol = 0:maxTime/100:maxTime;                   % Sparse time sampling
t_sol1 = time;                                   % Same time sampling as in simulation                    

y_sol = k_b / k_d_mean * (1 - exp(-k_d_mean * t_sol));
y_sol1 = k_b / k_d_mean * (1 - exp(-k_d_mean * t_sol1));

%% Graph ABK (stochastic) and deterministic trajectories
figure1 = figure('Name','Birth-Death: --> A -->','NumberTitle','off');
set(figure1,'Position',[1 1 500 450]);                      hold on;
figure1.PaperUnits = 'inches';
figure1.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.

p1 = plot(time,avgA,'b','Linewidth',2,'DisplayName','$\textrm{ABK } < \! N_A(t) \! >$');                 
p1_dev1 = plot(time,avgA+sdevA,'LineStyle','--','Color',[0.8 0.8 0.8]);
p1_dev0 = plot(time,avgA-sdevA,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p1s = plot(time,Sa(1,:),'m','MarkerSize',3,'DisplayName','Sample Trajectory');

p1d = plot(t_sol,y_sol,':b','DisplayName','$\textrm{DE } N_A(t)$','Linewidth',2);        

xlabel('$t \; \textrm{(sec)}$','Interpreter','LaTeX');           
ylabel('$N_A$','Interpreter','LaTeX');              hold off;   
% title([''],'FontName','Times New Roman','FontSize',11);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
axis([0 maxTime 0 15]);                                           
leg1 = legend([p1 p1d]);
% set(leg,'OuterPosition',[0.673 0.374 0.182 0.123]);
set(leg1,'Location','NorthWest');
set(leg1,'Interpreter','LaTeX','FontSize',10,...
    'EdgeColor',[0.95 0.95 0.95]);

% ***** Find Coefficient of Determination, R^2 *****
% This works for any curve fit, linear or nonlinear.
SST = sum((avgA - mean(avgA)).^2);        % Total sum of squares for simulation data (for A)
SSR = sum((avgA - y_sol1).^2);            % sum of square residuals (sim data vs DE predictions)

Rsq = 1 - SSR./SST                         % Definition of R^2
% R = sqrt(Rsq);                             % Correlation Coefficient R

% Determine Correlation Coefficient for <ABK> vs DE - Do pairwise comparison for each time point
% [CC, p_value] = corrcoef([y_sol1 avgA'])
% **************************************************

% Plot "instantaneous" mean k values (and their standard deviation) in an INSET PLOT
ax1_inset = axes('Parent',figure1,'Position',[0.56 0.19 0.25 0.25],'FontSize',7.5);        hold on;
p1_inset = plot(time(4:end),k_d_instavg(4:end),'Color',[0 0.6 0],...
    'LineWidth',2,'DisplayName','$\textrm{ABK } <k_d(t)>$','Parent',ax1_inset);
p1_inset_dev1 = plot(time(4:end),k_d_instavg(4:end)+k_d_inststd(4:end),...
    'LineStyle','--','Color',[0 1 0],'Parent',ax1_inset);
p1_inset_dev0 = plot(time(4:end),k_d_instavg(4:end)-k_d_inststd(4:end),...
    'LineStyle','--','Color',[0 1 0],'Parent',ax1_inset);
p1_inset_th = plot(time,k_d_mean*ones(1,size(time,2)),':r',...
    'LineWidth',1,'DisplayName','$\bar{k}_d$','Parent',ax1_inset);

axis([0 maxTime k_d_mean*0.6 k_d_mean*1.3]);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
set(gca,'XTick',0:100:500,'YTickLabel',{'0.006','0.008','0.010','0.012'});

xlabel('$t \; \textrm{(sec)}$','Interpreter','LaTeX');           
                                                                      hold off;
ylabel('$k_d \; \textrm{(sec} ^{-1} \textrm{)}$','Interpreter','LaTeX');                                        

% leg1_inset = legend(ax1_inset,[p1_inset p1_inset_th]);
% set(leg1_inset,'Position',[0.72 0.40 0.18 0.07],...
%     'Interpreter','LaTeX','FontSize',7,'EdgeColor',[0.95 0.95 0.95]);

% Create textbox
annotation(figure1,'textbox',...
    [0.80 0.28 0.20239121909187 0.045185184567063],'Color',[0 0.6 0],...
    'String','$\textrm{ABK } <\!k_d(t)\!>$',...
    'LineStyle','none','Interpreter','LaTeX','FontSize',9,'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
    [0.80 0.32 0.0564516129032256 0.0452962962962965],'Color',[1 0 0],...
    'String','$\bar{k}_d$',...
    'LineStyle','none','Interpreter','latex','FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
    [0.00771 0.93471 0.05010 0.06429],'String','a)',...
    'FontSize',12,'FontName','Times New Roman',...
    'FitBoxToText','off','EdgeColor','none');

%% Compare simulation standard deviation with CME prediction
% obtained from the Chemical Master Equation (see notes in this folder).
% Var = < n^2 > - < n >^2  = < n >           (True in this case; see derivation)
% SDev = sqrt(Var) = sqrt(<n>) 

% Plot SDev: ABK vs CME
figure2 = figure('Name','SDev Comparison','NumberTitle','off');
set(figure2,'Position',[501 1 500 450]);                          hold on;    
figure2.PaperUnits = 'inches';
figure2.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.

p2(1) = plot(time,sdevA,'LineWidth',2,'DisplayName','$\textrm{ABK}$');                          
p2(2) = plot(t_sol1,sqrt(avgA),'--','LineWidth',2,...
    'DisplayName','$\sqrt{ < \! N_A(t) \! > }$');
p2(3) = plot(t_sol,sqrt(y_sol),'-.','Color',[1 0.7 0.7],...
    'LineWidth',1,'DisplayName',...
    '$ \textrm{CME (using } \bar{k}_d \textrm{)} $');

% for a=1:size(time,2)
%     I(a) = integral(@(t) k_d_instavg(a) * t,0,time(a));
% end
% temp2 = sqrt(k_b ./ k_d_instavg .* (1 - exp(-I)));
temp = sqrt(k_b ./ k_d_instavg .* (1 - exp(- k_d_instavg .* time)));
plot(time,temp,'--k','DisplayName',...
    '$\sqrt{ \frac{k_b}{\bar{k}_d} \left( 1 - e^{-<k_d(t)> \, t} \right) }$');

xlabel('$t \; \textrm{(sec)}$','Interpreter','LaTeX');           
ylabel('$SDev \Big( <N_A(t)> \Big)$','Interpreter','LaTeX');                 
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
set(gca,'YTickLabel',{'0','0.5','1.0','1.5','2.0','2.5','3.0','3.5'});
axis([0 maxTime 0 3.5]);                                        hold off;

leg2 = legend(p2);
set(leg2,'Interpreter','LaTeX','EdgeColor',[0.95 0.95 0.95],'Position',[0.14 0.83 0.3 0.15])
% set(leg2,'FontName','Times New Roman','FontSize',9);

% *** Find Coefficient of Determination, R^2 ***
SST_s = sum((sdevA - mean(sdevA)).^2);          % Total sum of squares for simulation data (for A)
SSR_s = sum((sdevA - sqrt(y_sol1)).^2);         % sum of square residuals (sim data vs DE predictions)

Rsq_s = 1 - SSR_s ./ SST_s                      % Definition of R^2
% **********************************************

ax2_inset = axes('Parent',figure2,'Position',[0.6 0.19 0.25 0.25],'FontSize',7.5);        hold on;
p2_inset = plot(time,cvA,'k','DisplayName','ABK','LineWidth',2,'Parent',ax2_inset);               
pt2_inset(1) = plot(time,1./sqrt(avgA),':','Color',[0.2 0.75 0.5],...
    'DisplayName','Poisson','LineWidth',2,'Parent',ax2_inset);
% pt2_inset(2) = plot(t_sol,1./sqrt(y_sol),':','Color',[0.7 0.75 0.5],...
%     'DisplayName','Poisson, bar k','LineWidth',1,'Parent',ax2_inset);

axis([0 maxTime 0 1.0]);        
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
set(gca,'XTick',0:100:500);
set(gca,'YTick',0:0.2:1,'YTickLabel',{'0','0.2','0.4','0.6','0.8','1.0'});
leg2_inset = legend([p2_inset pt2_inset(1)]);
set(leg2_inset,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95]);
set(leg2_inset,'Position',[0.68 0.38 0.18 0.06]);
% set(leg2_inset,'Location','NorthEast');

xlabel('$t \; \textrm{(sec)}$','Interpreter','LaTeX');           
ylabel('$\eta$','Interpreter','LaTeX');                                                   hold off;                                    

% Create textbox
annotation(figure2,'textbox',...
    [0.00771 0.93471 0.05010 0.06429],'String','b)',...
    'FontSize',12,'FontName','Times New Roman',...
    'FitBoxToText','off','EdgeColor','none');

%% Plot instantaneous mean k values for all reps
% fig4 = figure('Name','Instantaneous k_d statistics','NumberTitle','off');                  
% set(fig4,'Position',[501 1 500 450]);                                     hold on;                  
% fig4.PaperUnits = 'inches';
% fig4.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.
% 
% for z=1:reps
%     p4(z) = plot(time,k_d_instavg_all(z,:));
% end
% p4_th = plot(time,k_d_mean*ones(1,size(time,2)),'r','LineWidth',5);
% 
% set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
% xlabel('t (sec)');           ylabel('k_d(t)');                           hold off;   

%% Finish
close(pb);
clear pb progress i j n t x temp tempA1 Av At p1* fig* leg* ax* SS*;
toc
