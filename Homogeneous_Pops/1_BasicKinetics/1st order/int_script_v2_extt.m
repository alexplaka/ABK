%  A --> B  ; Simulating 1st order kinetics. Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

clear; clc;             tic

% Assume constant time step increment (sec)
dt = 0.01;
totalTime = 10;                                 % Simulation time (sec)
t_steps = totalTime / dt;
time = 0:dt:totalTime-dt;

agents = 10;                                    % Initial number of A agents
lambda = 0.8;                   
k = log(2)/lambda;

reps = 1000;                                      % Number of times to repeat experiment

% k_calc = zeros(size(dt_array,2),n);
S = zeros(reps,t_steps);                        % To monitor population size
extt = zeros(reps,agents);                     % To monitor extinction time of agents

for n = 1:reps

    A = ones(1,agents);      
    S(n,1) = agents;
    P = zeros(1,t_steps);
    
    for t = 2:t_steps
        
        tempA = find(A==1);                     % find "alive" A agents
        
        for i = 1:size(tempA,2)
                
                P(t) = 1 - 2^(-dt/lambda);      % P from integrated rate law
%                     P(t) = k*dt;              % P from differential rate law

                if rand < P(t)                  % **** check probability condition ****
                    A(tempA(i)) = 0;            % agent dies
                    extt(n,tempA(i)) = time(t);        
                end         
                
        end                 % 'for i' loop
        
        S(n,t) = sum(A);
                
    end                     % 'for t' loop
    
%     k_calc(y,x) = FitCurve_exp_FixA(time',S',agents);
end                         % 'for n' loop

avg = mean(S);                          sdev = std(S);
y_theory = agents * exp(-k .* time);

%% Remove zeros from extt (zeros represent agents that didn't die within totalTime)
for z=1:agents
    extt_zeros = find(extt(:,z)==0);
    for a = 1:size(extt_zeros,1)
        extt(extt_zeros(a),z) = NaN;
    end
end
extt = reshape(extt,[],1);                      % Reshape into column vector

%% Plot results
figure;
p1 = plot(time,avg,'b');           hold on;
p2 = plot(time,avg+sdev,'k','LineStyle','--','Color',[0.8 0.8 0.8]);
p3 = plot(time,avg-sdev,'k','LineStyle','--','Color',[0.8 0.8 0.8]);
p4 = plot(time,y_theory,':g');

title(['N_{A,i} = ' num2str(agents)],'FontName','Times New Roman','FontSize',12);
axis([0 totalTime 0 agents]);
set(gca,'XMinorTick','on','Box','off');
xlabel('t (sec)');                      ylabel('N_A');

leg0 = legend([p1 p4],'<N_A(t)>_{sim}','N_A(t) = 10\cdote^{-k\cdott}');
set(leg0,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');       

% hgsave(['o1_TCdev_A=' num2str(agents) '.fig']);
%% Plot Extinction time distribution 
bins = 20;

figure;
% hist(extt);                                   % histogram for all agents
% h = findobj(gca,'Type','patch');              % use set(h(1),...) commands below to customize

[num1, c1] = hist(extt,bins);                   % obtain histogram data
Num1 = num1 / (reps*agents);                    % Normalize

b = bar(c1',Num1');                                             hold on;
set(b(1),'FaceColor','flat','EdgeColor',[0.70 0.78 1]);  

exp_PDF1 = 1 / (bins/totalTime) * k * exp(- k .* time);
% PDFs need to be adjusted for the number of bins in the histogram
% in order to obtain proper fits. Adjustment: divide by (bins/totalTime)

plot(time,exp_PDF1,'-','Color','g','LineWidth',2);

% *** Use curve-fitting to verify Exponential Dist PDF. ***
% [est1,y_fit1] = FitCurve_exp_Fixk(c1',Num1',k);       % fit for fixed k
% [est1,y_fit1] = FitCurve_exp(c1',Num1');
% * Graph Fits *
% plot(c1,y_fit1,'-','Color',[0.15 0.15 1],'LineWidth',1);
% *********************************************************

set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
% set(gca,'XTick',[0 1 2 3 4 5 6 7 8 9 10]);
xlabel('t (sec) ');             ylabel('Fraction of A\rightarrowX Transitions');

% leg1 = legend(['Subspecies A1, k_{A1} = ' num2str(k_sub(1)) ' sec^{-1}']);
% set(leg1,'FontName','Times New Roman','FontSize',9,...
%     'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');    

% Coefficients of Determination R^2 for histogram fits to exponential dist PDF.
exp_PDF1d = 1 / (bins/totalTime) * k * exp(- k * c1);       % Discretize exp PDF

disp(['R^2 = ' num2str(CoefDet(Num1,exp_PDF1d))]);

% ** Note: the distribution is the same no matter what the initial population is. **
%% Finish
toc
clear a i t x z p* leg* P A tempA;
save(['o1_A=' num2str(agents) '_extt.mat']);
