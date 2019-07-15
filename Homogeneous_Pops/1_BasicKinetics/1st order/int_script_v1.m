%  A --> B  ; Simulating 1st order kinetics. Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

clear; clc;             tic

% Assume constant time step increment (sec)
dt = 0.01;
totalTime = 10;         % Simulation time (sec)
t_steps = totalTime / dt;

agents = 10;

lambda = 0.8;                   k = log(2)/lambda;

reps = 1000;                            % Number of times to repeat experiment

% k_calc = zeros(size(dt_array,2),n);
S = zeros(reps,t_steps);

P = 1 - 2^(-dt/lambda);      % P from integrated rate law (constant)
%  P = k*dt;                    % P from differential rate law (constant)

for n = 1:reps
    fprintf('.');
    
    A = ones(1,agents);      
    S(n,1) = agents;
    time = zeros(1,t_steps);
    
    for t = 2:t_steps
        
        tempA = find(A==1);                     % find "alive" A agents
                
        for i = 1:size(tempA,2)
                if rand < P                     % ****check probability condition
                    A(tempA(i)) = 0;            % agent dies
                end         % 'if rand' statement
        end                 % 'for i' loop
        S(n,t) = sum(A);
    end                     % 'for t' loop
    
%     k_calc(y,x) = FitCurve_exp_FixA(time',S',agents);
end                         % 'for n' loop

clear t_steps n t i tempA x A;

avg = mean(S);                      sdev = std(S);              cv = sdev ./ avg;
time = 0:dt:totalTime-dt;

%% Theoretical time trajectory
t_theory = 0:dt:totalTime;                   y_theory = agents*exp(-k.*t_theory);

%% Plot time trajectories
fig1 = figure;
set(fig1,'Position',[1 1 500 450]);                          hold on;    
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.

p1 = plot(time,avg,'b');           hold on;
p2 = plot(time,avg+sdev,'k','LineStyle','--','Color',[0.8 0.8 0.8]);
p3 = plot(time,avg-sdev,'k','LineStyle','--','Color',[0.8 0.8 0.8]);
p4 = plot(t_theory,y_theory,':g');
set(gca,'XMinorTick','on','Box','off');
title(['N_{A,i} = ' num2str(agents)],'FontName','Times New Roman','FontSize',12);
axis([0 10 0 agents]);
xlabel('t (sec)');                      ylabel('N_A');

leg1 = legend([p1 p4],'ABK <N_A(t)>',['N_A(t) = ' num2str(agents) ' e^{-kt}']);
set(leg1,'Location','SouthEast');
set(leg1,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95]);

% Create textbox
annotation(fig1,'textbox',...
    [0.00771 0.93471 0.05010 0.06429],'String','b)',...
    'FontSize',12,'FontName','Times New Roman',...
    'FitBoxToText','off','EdgeColor','none');

% hgsave(['o1_TCdev_A=' num2str(agents) '.fig']);
%% Plot coefficient of variation
fig2 = figure;                                                  hold on;
% set(fig2,'Position',[1 1 500 450]);                          
fig2.PaperUnits = 'inches';
fig2.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.

p = plot(time,cv,'DisplayName','ABK');
pt = plot(t_theory,1./sqrt(y_theory),'DisplayName','Poisson');

% set([p(2) p(3) pt(2) pt(3)],'Visible','off');       % Don't show for ...
% axis([0 time(end) 0 1]);        
axis tight

leg2 = legend([p pt(1)]);
set(leg2,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95],'Location','Best');

xlabel('t (sec)');         
ylabel('$\eta$','Interpreter','LaTeX');
hold off;

%% Finish
clear fig* leg* p* *theory;
mkdir(['reps=' num2str(reps)]);
cd(['reps=' num2str(reps)]);
save(['o1_A=' num2str(agents) '.mat']);
cd ..
toc