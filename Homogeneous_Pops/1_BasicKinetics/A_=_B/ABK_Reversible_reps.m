% A <--> B  
% Simulating 1st order kinetics for reversible chemical reaction.
% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;            tic;
rng(0);
pb = waitbar(0,'0');
%% Declare variables and functions
global ko_f ko_r;

maxTime = 25;                        % Maximum Simulation time (sec)
dt = 1/100;                             % Constant (fixed) time step increment (sec)
t_steps = maxTime / dt;

agents = 20;          	disp(['Agents = ' num2str(agents)]);

ko_f = 0.04;                               % basal forward rate
ko_r = 0.20;                               % basal reverse rate

reps = 500;
Sa = zeros(reps,t_steps);                 Sb = zeros(reps,t_steps);         
%% ABK Simulation 

for n=1:reps
    progress = n/reps;
    waitbar(progress,pb,sprintf('%.0f%% %',progress*100));

    % Initialize - Preallocate memory for variables; ** Initial conditions **
    P = zeros(1,t_steps);               % For storing probability value at each time step
    At = zeros(1,t_steps);              % Sum of A agents in each time step	
    Bt = zeros(1,t_steps);              % Sum of B agents in each time step

    % ******** Initial conditions - Number of A, B Agents ********
    At(1) = floor(agents/2);            Bt(1) = agents - At(1);					
    % ************************************************************

    tempA = zeros(1,agents);            % tempB = zeros(1,agents);              
    % Put "1" where agents are "alive", then randomize the array
    for c=1:At(1),                      tempA(c)=1;             end
    tempA = RandArray(tempA);           % Randomize array
    Av = [tempA ; tempA];               % Initialize vector storing state of A agents     
    Bv = ~Av;                           % Initialize vector storing state of B agents 
    % Notes on Av, Bv:
    % - Markov process, so only previous and current time steps needed --> 2 rows:          
    % - Av and Bv are complementary (A+B=agents)
    clear c tempA tempB;

    time = zeros(1,t_steps);
    
    for t=2:t_steps
    %     dt = 1 / (ko*Sa(t-1));            % Variable time step increment ; OLD
    %     dt = exprnd(1/(k*S(t-1)));   % Exponentially-distributed time step increment ; OLD

        for i = 1:agents
            if Av(1,i) == 1                     % if A agent is alive
    %             P(t-1) = 1 - exp(-ko_f * dt);      % P from integrated rate law
                P(t-1) = ko_f * dt;                % P from differential rate law
                if rand < P(t-1)                  % **check probability condition**
                    Av(2,i) = 0;                 % A agent dies
                    Bv(2,i) = 1;                 % A is converted to B
                end
            elseif Bv(1,i) == 1					% if B agent is alive
    %        	  P(t-1) = 1 - exp(-ko_r * dt);      % P from integrated rate law
                P(t-1) = ko_r * dt;                % P from differential rate law
                if rand < P(t-1)                  % **check probability condition**
                    Bv(2,i) = 0;                 % B agent dies
                    Av(2,i) = 1;                 % B is converted to A
                end
            end
        end
        At(t) = sum(Av(2,:));                  Bt(t) = sum(Bv(2,:));
        Av(1,:) = Av(2,:);                     Bv(1,:) = Bv(2,:); 
        time(t) = time(t-1) + dt;
    end
    Sa(n,:) = At;                       Sb(n,:) = Bt;
end         % end "for n" loop

avgA = mean(Sa);                     sdevA = std(Sa);
avgB = mean(Sb);                     sdevB = std(Sb);

A_ss = ko_r / (ko_f + ko_r) * agents;
disp(['A_ss = ' num2str(A_ss)]);
%% Solve differential equation (was waiting for variable 'finaltime')
[t_sol, y_sol] = ode45(@o1_rev_dif,time,[At(1) ; Bt(1)]);     % Same time sampling as in simulation
% Or, uncomment the following to try a different time sampling
% [t_sol, y_sol] = ode45(@o1_rev_dif,0:maxTime/100:maxTime,[At(1) ; Bt(1)]);

disp(['Diff eq terminal A value   = ' num2str(y_sol(end,1))]);

%% Find Coefficient of Determination, R^2
% This works for any curve fit, linear or nonlinear.
SST(1) = sum((avgA' - mean(avgA)).^2);     % Total sum of squares for simulation data (for A)
SSR(1) = sum((avgA' - y_sol(:,1)).^2);     % sum of square residuals (sim data vs DE predictions)

SST(2) = sum((avgB' - mean(avgB)).^2);     % Total sum of squares for simulation data (for B)
SSR(2) = sum((avgB' - y_sol(:,2)).^2);     % sum of square residuals (sim data vs DE predictions)

Rsq = 1 - SSR./SST                         % Definition of R^2
R = sqrt(Rsq);                             % Correlation Coefficient R

%% Graph ABM (stochastic) and deterministic results
figure1 = figure('Name','A <--> X','NumberTitle','off');
set(figure1,'Position',[1 1 500 450]);                  hold on;
p1 = plot(time,avgA,'b','MarkerSize',3,'DisplayName','<N_A(t)>_{sim}');                 
% p1_dev1 = plot(time,avgA+sdevA,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p1_dev0 = plot(time,avgA-sdevA,'LineStyle','--','Color',[0.8 0.8 0.8]);
p1s = plot(time,Sa(end-2,:),'m','MarkerSize',3,'DisplayName','Sample Trajectory');

p2 = plot(time,avgB,'r','MarkerSize',3,'DisplayName','<N_B(t)>_{sim}');                 
p2_dev1 = plot(time,avgB+sdevB,'LineStyle','--','Color',[0.8 0.8 0.8]);
p2_dev0 = plot(time,avgB-sdevB,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p2s = plot(time,Sb(end,:),'m','MarkerSize',3,'DisplayName','Sample Trajectory');

p1d = plot(t_sol,y_sol(:,1),':b','DisplayName','DE');        
p2d = plot(t_sol,y_sol(:,2),':r');

xlabel('t (sec)');           ylabel('N(t)');             hold off;   
axis([0 maxTime 0 agents]);                                           
% title(['k_r/k_f = ' num2str(ko_r/ko_f) '/1'],...
%     'FontName','Times New Roman','FontSize',11);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
leg = legend([p1 p2]);
% set(leg,'OuterPosition',[0.673 0.374 0.182 0.123]);
set(leg,'Location','NorthWest');
set(leg,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95]);

%% Finish
close(pb);
clear pb progress r w i j temp Av Bv fig* p* leg;
toc

%% Notes
% - Works?  Yes, perfectly well.