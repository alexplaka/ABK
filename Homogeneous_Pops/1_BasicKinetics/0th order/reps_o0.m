%  --> X                    
% Simulating 0th order kinetics.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

clear; clc;                                 tic

k = 1;      % Zeroth order MICROSCOPIC kinetic constant [units: N_x / sec];

% agents = 1000;          % Max number of A molecules that can be produced 

t_max = 100;                            % in seconds
dt = 1/100;                             % Fixed time step increment
steps = t_max / dt;
time = zeros(1,steps);

Xi = 0;                     % Initial population size of A

reps = 10000;
Sx = zeros(reps,steps);

for n=1:reps
    
    X = zeros(1,steps);
    
    for t=2:steps               

        P = k * dt;               % Note: Does NOT depend on A
                 
        if rand < P       
            X(t) = X(t-1) + 1;
        else
            X(t) = X(t-1);
        end

        time(t) = time(t-1) + dt;
    end
    Sx(n,:) = X;
end

Sx = Sx + Xi; 

avg = mean(Sx);                     sdev = std(Sx);

%% ODE and CME solutions for 0th order kinetics.
DE_sol = k * time + Xi;                 % ODE solution
% Note: above ODE solution is the same as the CME predicts for the average trajectory.
CME_sdev = sqrt(DE_sol);                % CME-predicted standard deviation

% Find Coefficient of Determination, R^2
Rsq_avg = CoefDet(avg, DE_sol)
Rsq_sdev = CoefDet(sdev, CME_sdev)

%% Plot time course
fig1 = figure('Name','0th Order Rx Time Course','NumberTitle','off');         
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.
fig1.Position = [1 1 500 450];

hold on;

p1 = plot(time,avg,'b','MarkerSize',3,...
    'DisplayName','$\textrm{ABK } < \! N_X(t) \! >$');             
p1_dev1 = plot(time,avg+sdev,'LineStyle','--','Color',[0.8 0.8 0.8]);
p1_dev0 = plot(time,avg-sdev,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p1_de = plot(t_sol,y_sol(:,1),':r','DisplayName','$N_X(t)=k_0 \, t$');
p1_de = plot(time,DE_sol,':r','MarkerSize',3,...
    'DisplayName','$N_X(t)=k_0 \, t$');

trial = 11;                                 % Should be <= reps
p1_s = plot(time,Sx(trial,:),'g','MarkerSize',3,...
    'DisplayName','$\textrm{Sample Trajectory}$');

axis tight;  
axis([0 round(time(end)) 0 round(max(avg+sdev+1))]);  
hold off;
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('$t \; \textrm{(sec)}$','Interpreter','LaTeX','FontSize',11);                 
ylabel('$N_X(t)$','Interpreter','LaTeX','FontSize',11); 

leg1 = legend([p1 p1_de p1_s]);
set(leg1,'Location','NorthWest','Interpreter','LaTeX',...
    'EdgeColor',[0.95 0.95 0.95],'FontSize',10);

%% Plot standard deviation 
fig2 = figure('Name','0th Order Rx Standard Deviation','NumberTitle','off');         
fig2.PaperUnits = 'inches';
fig2.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.
fig2.Position = [1 1 500 450];      

hold on;
p2_ABK = plot(time,sdev,'DisplayName','ABK','LineWidth',2);
p2_CME = plot(time,CME_sdev,'--','DisplayName','CME','LineWidth',2);
hold off;

axis([0 round(time(end)) 0 round(sqrt(k*time(end))+1)]);                
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('$t \; \textrm{(sec)}$','Interpreter','LaTeX','FontSize',11);                 
ylabel('$SDev \Big( < N_X(t) > \Big) $','Interpreter','LaTeX','FontSize',11); 

leg2 = legend([p2_ABK p2_CME]);
set(leg2,'Location','NorthWest','Interpreter','LaTeX',...
    'EdgeColor',[0.95 0.95 0.95],'FontSize',10);


%% Finish
clear n t p* A fig* leg* SS* R;
toc
% Result: It works! 
