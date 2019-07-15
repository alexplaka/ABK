%  A --> B  ; Simulating 1st order kinetics. Considering individual agents.
% "extt" is short for agent-EXTinction-Time.
% This script is checking that the distribution of "death" events is a
% decaying exponential with rate constant k.
% Also see similar treatment for 1st order process in heterogeneous pops.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

clear; clc;             tic

% Assume constant time step increment (sec)
dt = 0.01;
totalTime = 10;         % Simulation time (sec)
t_steps = totalTime / dt;
agents = 10;
k = 0.7;
reps = 500;              % Number of times to repeat experiment
% k_calc = zeros(size(dt_array,2),reps);
S = zeros(reps,t_steps);
extt = [];                  % For recording the time of "extinction/death" events

for x = 1:reps

    A = ones(t_steps,agents);      
    S(1) = agents;
    time = zeros(1,t_steps);
    
    P = 1 - exp(-k*dt);                 % P from integrated rate law
    % P = k*dt;                           % P from differential rate law
    
    for t = 2:t_steps
        time(t) = time(t-1) + dt;
        for i = 1:size(A,2)
            if A(t,i) == 1                 % if agent is still alive
                if rand < P                  % ****check probability condition
                    A(t,i) = 0;            % agent dies
                    extt = [extt ; time(t)];
                    for w = t+1:t_steps
                        A(w,i) = 0;        % mark agent dead for rest of time
                    end     % 'for w' loop
                end         % 'if rand' statement
            end             % 'if A' statement
        end                 % 'for i' loop
        S(x,t) = sum(A(t,:));
    end                     % 'for t' loop
    clear t w i;
%     k_calc(y,x) = FitCurve_exp_FixA(time',S',agents);
end                         % 'for x' loop
clear t_steps x P A;

avg = mean(S);                          sdev = std(S);
t = 0:0.01:totalTime;                   y_theory = agents*exp(-k.*t);

%% Plot results
figure;
p1 = plot(time,avg,'b');           hold on;
p2 = plot(time,avg+sdev,'k','LineStyle','--','Color',[0.8 0.8 0.8]);
p3 = plot(time,avg-sdev,'k','LineStyle','--','Color',[0.8 0.8 0.8]);
p4 = plot(t,y_theory,':g');
set(gca,'XMinorTick','on','Box','off');
legend([p1 p4],'<N_A>_{sim}','N_A = 10\cdote^{-k\cdot\Deltat}')
title(['N_{A,i} = ' num2str(agents)],'FontName','Times New Roman','FontSize',12);
axis([0 10 0 agents]);
xlabel('t (sec)');                      ylabel('N_A');
% hgsave(['o1_TCdev_A=' num2str(agents) '.fig']);

%% Histogram
figure;
hist(extt);                             % histogram of death events
bins = 20;
[num c] = hist(extt,bins);              % histogram data
Num = num / (reps*agents);              % Normalize
[est , y_fit] = FitCurve_exp_Fixk(c',Num',k);           % Fixed k
bar(c,Num);
hold on;
plot(c,y_fit,'-','Color',[0.5 0.1 1],'LineWidth',1);
set(gca,'XTick',[0 1 2 3 4 5 6 7 8 9 10],...
    'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('t (sec) ');             ylabel('Fraction of "Death" Events');

%% Finish
% save(['o1_A=' num2str(agents) '.mat']);
toc