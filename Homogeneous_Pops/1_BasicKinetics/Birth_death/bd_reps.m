%  --> A -->  
%  k_b   k_d
% Simulating birth-death process for A: 
% k_b: birth, 0th order process; k_d: death, 1st order process. 
% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;            tic;
rng(1);
pb = waitbar(0,'0');
%% Declare variables and functions
global k_b k_d;

k_b = 0.1;                    % Microscopic 0th order birth rate (units: 1/sec)
k_d = 0.01;                   % "Death" rate; 1st order (units: 1/sec)

A_ss = k_b / k_d;             % Steady-state value
disp(['Predicted A_ss = ' num2str(A_ss)]);

agents = floor(3 * A_ss);
% agents = 100;

% ******** Initial condition - Number of A Agents ********
Ao = 0;                       % Initial population size of A

maxTime = 500;                % Maximum Simulation time (sec)  ********************
dt = 1/100;                   % Constant (fixed) time step increment (sec)
t_steps = maxTime / dt;

% ***** Probability condition can be specified here; they are constant ******
% P_b = k_b * dt;               % Probability of "birth" process (0th order), P_dif
P_b = 1 - exp(-k_b*dt);       % Probability of "birth" process (0th order), P_ber

% P_d = k_d * dt;               % Probability of "death" process (1st order), P_dif
P_d = 1 - exp(-k_d * dt);     % Probability of "death" process (1st order), P_ber
% *********************************************************************************

reps = 500;  % change this to 5000 to obtain figure 2.15b,c in the report.
Sa = zeros(reps,t_steps);    
%% ABK Simulation 

for n=1:reps
    
    progress = n/reps;
    waitbar(progress,pb,sprintf('%.0f%% %',progress*100));

    % Initialize - Preallocate memory for variables; ** Initial conditions **
    At = zeros(1,t_steps);              % Sum of A agents in each time step	

    At(1) = Ao;          				
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

        if rand < P_b                           % "Birth", 0th order reaction
            temp = find(Av(1,:)==0);            % Randomly choose agent which becomes "alive" 
            Av(2,temp(ceil(rand * size(temp,2)))) = 1;
        end

        tempA1 = find(Av(1,:)==1);
        for i = 1:size(tempA1,2)
            if rand < P_d                       % "Death", 1st order reaction
                Av(2,tempA1(i)) = 0;            % A agent dies
            end
        end

        At(t) = sum(Av(2,:));                  
        Av(1,:) = Av(2,:);                     
        time(t) = time(t-1) + dt;
    end
    
    Sa(n,:) = At;       
    
end         % end "for n" loop

% See if extinction (N_A(t)=0) occurred in any of the trajectories
for x=1:size(Sa,1)
    if Sa(x,end)==0
        disp(['rep = ' num2str(x) ' ; Extinction!']);
    end
end

%% Calculate Average and SDev of time trajectory   (* Start from here when loading data *)
avgA = mean(Sa);                     sdevA = std(Sa);               
cvA = sdevA ./ avgA;                                % Coefficient of variation

%% Compute deterministic trajectory (using analytic solution)
t_sol = 0:maxTime/100:maxTime;                   % Sparse time sampling
t_sol1 = time;                                   % Same time sampling as in simulation                    

y_sol = k_b / k_d * (1 - exp(-k_d*t_sol));
y_sol1 = k_b / k_d * (1 - exp(-k_d*t_sol1));

%% Find Coefficient of Determination, R^2
% This works for any curve fit, linear or nonlinear.
SST = sum((avgA - mean(avgA)).^2);        % Total sum of squares for simulation data (for A)
SSR = sum((avgA - y_sol1).^2);            % sum of square residuals (sim data vs DE predictions)

Rsq = 1 - SSR./SST                         % Definition of R^2
% R = sqrt(Rsq);                             % Correlation Coefficient R

% Determine Correlation Coefficient for <ABK> vs DE - Do pairwise comparison for each time point
% [CC, p_value] = corrcoef([y_sol1 avgA'])

%% Graph ABK (stochastic) and deterministic trajectories
figure1 = figure('Name','Birth-Death: --> A -->','NumberTitle','off');
set(figure1,'Position',[1 1 500 450]);                      hold on;
figure1.PaperUnits = 'inches';
figure1.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.

p1 = plot(time,avgA,'b','MarkerSize',3,'DisplayName','<N_A(t)>_{sim}');                 
p1_dev1 = plot(time,avgA+sdevA,'LineStyle','--','Color',[0.8 0.8 0.8]);
p1_dev0 = plot(time,avgA-sdevA,'LineStyle','--','Color',[0.8 0.8 0.8]);
p1s = plot(time,Sa(end-11,:),'m','MarkerSize',3,'DisplayName','Sample Trajectory');

p1d = plot(t_sol,y_sol,':b','DisplayName','DE N_A(t)');        

xlabel('t (sec)');           ylabel('N_A');              hold off;   
% title([''],'FontName','Times New Roman','FontSize',11);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
axis([0 maxTime 0 15]);                                           
leg1 = legend([p1 p1d p1s]);
% set(leg,'OuterPosition',[0.673 0.374 0.182 0.123]);
set(leg1,'Location','SouthEast');
set(leg1,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95]);

% Create textbox
annotation(figure1,'textbox',...
    [0.00771 0.93471 0.05010 0.06429],'String','a)',...
    'FontSize',12,'FontName','Times New Roman',...
    'FitBoxToText','off','EdgeColor','none');

%% Compare simulation standard deviation with deterministic expectation...
% obtained from the Chemical Master Equation (see notes in this folder).
% Var = < n^2 > - < n >^2  = < n >           (True in this case; see derivation)
% SDev = sqrt(Var) = sqrt(<n>) 

% Plot SDev: ABK vs CME
figure2 = figure('Name','SDev Comparison','NumberTitle','off');
set(figure2,'Position',[1 1 500 450]);                          hold on;    
figure2.PaperUnits = 'inches';
figure2.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.

plot(time,sdevA,'LineWidth',2);                          
plot(t_sol,sqrt(y_sol),'--','LineWidth',2);
xlabel('t (sec)');           
ylabel('$SDev \Big( <N_A(t)>_{sim} \Big)$','Interpreter','LaTeX');                 
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
set(gca,'YTickLabel',{'0','0.5','1.0','1.5','2.0','2.5','3.0','3.5'});
axis([0 maxTime 0 3.5]);                                        hold off;

leg2 = legend('ABK','CME');
set(leg2,'Location','SouthEast');
set(leg2,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95]);

% Create textbox
annotation(figure2,'textbox',...
    [0.00771 0.93471 0.05010 0.06429],'String','b)',...
    'FontSize',12,'FontName','Times New Roman',...
    'FitBoxToText','off','EdgeColor','none');

% Find Coefficient of Determination, R^2
SST_s = sum((sdevA - mean(sdevA)).^2);          % Total sum of squares for simulation data (for A)
SSR_s = sum((sdevA - sqrt(y_sol1)).^2);         % sum of square residuals (sim data vs DE predictions)

Rsq_s = 1 - SSR_s ./ SST_s                      % Definition of R^2

%% Plot coefficient of variation
fig3 = figure('Name','Coefficient of Variation','NumberTitle','off');                  
set(fig3,'Position',[1 1 500 450]);                                     hold on;                  
fig3.PaperUnits = 'inches';
fig3.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.

p = plot(time,cvA,'k','DisplayName','ABK','LineWidth',2);               
pt = plot(t_sol,1./sqrt(y_sol),':','Color',[0.2 0.75 0.5],'DisplayName','Poisson','LineWidth',2);

% set([p(2) p(3) pt(2) pt(3)],'Visible','off');       % Don't show for ...
axis([0 maxTime 0 1.0]);        
% axis tight;
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');

leg3 = legend([p pt(1)]);
set(leg3,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95],'Location','Best');

xlabel('t (sec)');         
ylabel('$\eta$','Interpreter','LaTeX');                                 hold off;                                    

%% Finish
close(pb);
clear pb progress i j n t x temp tempA1 Av At p1* fig* leg* SS*;
toc
