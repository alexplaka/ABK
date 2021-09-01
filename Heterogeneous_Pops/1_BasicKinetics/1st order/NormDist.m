%  A --> B  ; Simulating 1st order kinetics. Considering individual agents.
% *** Heterogeneous population of A ***
% Here, we assume normally distributed k values of the A population.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear; clc;             tic
rng(0);

% Assume constant time step increment (sec)
dt = 0.01;
totalTime = 10;         % Simulation time (sec)
t_steps = totalTime / dt;

reps = 500;              % Number of times to repeat experiment
% k_calc = zeros(size(dt_array,2),n);
S = zeros(reps,t_steps);

agents = 20;            % Population size (should be divisible by # of subspecies)

extt = zeros(reps, agents);     % To monitor extinction time of agents

% ****** Set population's k values ************
% Normally-distributed k values across population 
% with mean = ko and standard deviation = sigma.
ko = 0.5;                           sigma = 0.1;
repeat = 1;
while repeat == 1
    k = ko + sigma * randn(1,agents);  
    negk = find(k<=0);
    if isempty(negk)==1             % Ensure there are no negative k values
        repeat = 0;
    end
end
% # of subspecies = # of agents
% *********************************************

% Time-independent probability for a given agent
P = 1 - exp(-k .* dt);                  % P from integrated rate law
% P = k .* dt;                            % P from differential rate law

for n = 1:reps

    A = ones(2,agents);     % Only 2 rows because of Markov property      
    S(n,1) = agents;
    time = zeros(1,t_steps);
    
    for t = 1:t_steps
        tempAa = find(A(1,:)==1);       % find "alive" A agents
        
        for i = 1:size(tempAa,2)
            if rand < P(tempAa(i))              % ** check probability condition **
                A(2,tempAa(i)) = 0;            % agent dies
                extt(n,tempAa(i)) = time(t)+dt;
            end             % 'if rand' statement
        end                 % 'for i' loop
        
        S(n,t+1) = sum(A(2,:));
        
        A(1,:) = A(2,:);
        time(t+1) = time(t) + dt;
    end                     % 'for t' loop
end                         % 'for n' loop

avg = mean(S);                          sdev = std(S);

% Calculate theoretical curves
t = 0:dt:totalTime;                     y_theory = agents*exp(- ko * t);

% Remove zeros from extt (zeros represent agents that didn't die with totalTime)
for z=1:agents
    extt_zeros = find(extt(:,z)==0);
    for a = 1:size(extt_zeros,1)
        extt(extt_zeros(a),z) = NaN;
    end
end
%% Plot results
figure;
p1 = plot(time,avg,'b','LineWidth',1);           hold on;
p2 = plot(time,avg+sdev,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);
p3 = plot(time,avg-sdev,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);
p4 = plot(t,y_theory,':g','LineWidth',2);

set(gca,'XMinorTick','on','Box','off');
title(['N_{A,i} = ' num2str(agents)],'FontName','Times New Roman','FontSize',12);
axis([0 10 0 agents]);
xlabel('t (sec)');                      ylabel('N_A');

leg0 = legend([p1 p4],'<N_A(t)>_{sim}','N_A(t) = 10\cdote^{-k\prime\cdot\Deltat}');
set(leg0,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');       

% hgsave(['o1_TCdev_A=' num2str(agents) '.fig']);
%% Plot Extinction time distribution for each agent
figure;
hist(extt(:,:));           % Choose agents to compare in column entry
set(gca,'XTick',[0 1 2 3 4 5 6 7 8 9 10],...
    'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('t (sec) ');             ylabel('Number of "Death" Events');
%% Finish
clear a z n i;
% save(['o1_A=' num2str(agents) '.mat']);
toc