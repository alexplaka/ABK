%  A + B --> C
% Simulating 2nd order kinetics using ABK.
% The experiment is repeated in order to obtain average trajectories

% Heterogeneous Population.
% Here, I evaluate each possible interaction wrt A.
% (ie, I draw a random number for each possible A-B pair).

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear; clc;     tic
rng(0);

agents = 10;

maxTime = 200;
dt = 1/100;                            % Fixed time step increment
t_steps = ceil(maxTime / dt);

reps = 100;                              % Repeat experiment this # of times

% Initial number of reactant molecules; assume B is limiting.
Ao = agents;                    Bo = ceil(0.7*agents);                
Co = 0;                         Cmax = Bo;

% **************** 2 Distinct Subspecies *********************
% ********* Set up k-matrix wrt A (size: Ao x Bo) ************
subspecies = 2;
k1 = 0.001*ones(Bo,Ao/2);                   k2 = 0.003*ones(Bo,Ao/2);
k = [k1, k2]';
ko = mean(mean(k));
% ************************************************************

% **************** # Subspecies = Agents *********************
% **************** Normally-Distributed k ********************
% ********* Set up k-matrix wrt A (size: Ao x Bo) ************
% ko = 0.001;                     % Mean k value
% sigma = 0.0001;                 % standard deviation of k
% 
% repeat = 1;
% while repeat == 1
%     k = ko + sigma*randn(Ao,Bo);    % Normally distributed with mean=ko and std=sigma
%     repeat = 0;
%     negk = find(k<=0);
%     if isempty(negk)==0             % Ensure there are no negative k values
%         repeat = 1;
%     end
% end
% ************************************************************

extt = zeros(reps, Ao);     % To monitor extinction time of A agents

%% Preallocate memory; Initialize time-dependent sum of molecule numbers
At = zeros(reps,t_steps+1);        Bt = zeros(reps,t_steps+1);        Ct = zeros(reps,t_steps+1);

for n=1:reps

    fprintf(1,'.');
    
    % Arrays for tracking reactant and product agent status
    A = ones(1,Ao);                 B = ones(1,Bo);                 C = zeros(1,Cmax); 
    At(n,1) = sum(A);               Bt(n,1) = sum(B);               Ct(n,1) = sum(C);
    time = zeros(1,t_steps);
    
    for t=1:t_steps

        for i = 1:size(A,2)
            if A(i) == 1                    % if this A agent is alive
                for j = 1:size(B,2)         % then see if it reacts with a B agent
                    if B(j) == 1            % if this B agent is alive
                        P = 1 - exp(- k(i,j) * dt);
                        if rand < P
                            A(i) = 0;       
                            B(j) = 0;
                            C(j) = 1;  
                            
                            extt(n,i) = time(t)+dt;
                        end
                    end
                end
            end
        end 

        
    At(n,t+1) = sum(A);         Bt(n,t+1) = sum(B);         Ct(n,t+1) = sum(C);
    time(t+1) = time(t) + dt;
    
    end

end                 % end "for reps" loop

% Population statistics of time trajectories
avgA = mean(At);            avgB = mean(Bt);            avgC = mean(Ct);
sdevA = std(At);            sdevB = std(Bt);            sdevC = std(Ct);

% Solve ODE for 2nd order kinetics
[t_sol,y_sol] = ode45(@(t,y) o2_dif(t,y,ko),[0:maxTime/100:maxTime],[Ao ; Bo ; Co]);

%% ** Group extinction times by subspecies ** - Works for subspecies = 2 simulation 
% Remove zeros from extt (zeros represent agents that didn't die with totalTime)
f = size(extt,2);
for z=1:f
    extt_zeros = find(extt(:,z)==0);
    for a = 1:size(extt_zeros,1)
        extt(extt_zeros(a),z) = NaN;
    end
end

extt2 = zeros(reps*f/subspecies,subspecies);
s = 1;              w = 1;
for j=1:f
    extt2((w-1)*reps+1:w*reps,s) = extt(:,j);
    
    if j < f && k(j,1) ~= k(j+1,1)
        s = s + 1;
        w = 1;
    else
        w = w + 1;
    end  
end

%% Plot
figure('Name','2nd Order Rx Time course','NumberTitle','off'); 
plot(time,avgA,time,avgB,time,avgC);                                    hold on;
p1_devA1 = plot(time,avgA+sdevA,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);
p1_devA0 = plot(time,avgA-sdevA,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);
p1_devB1 = plot(time,avgB+sdevB,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);
p1_devB0 = plot(time,avgB-sdevB,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);
p1_devC1 = plot(time,avgC+sdevC,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);
p1_devC0 = plot(time,avgC-sdevC,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);

% Plot DE Time Courses
plot(t_sol,y_sol(:,1),':c','LineWidth',2);               % scatter(t_sol,y_sol(:,1),'.b');                  
plot(t_sol,y_sol(:,2),':g','LineWidth',2);               % scatter(t_sol,y_sol(:,2),'.g');
plot(t_sol,y_sol(:,3),':k','LineWidth',2);               % scatter(t_sol,y_sol(:,3),'.r');

title(['N_{A,i} = ' num2str(Ao) ' , N_{B,i} = ' num2str(Bo) ' , N_{X,i} = ' num2str(Co)],...
    'FontName','Times New Roman','FontSize',12);
xlabel('t (sec)');                  ylabel('Population Size');          hold off;
axis([0 maxTime 0 agents]);
set(gca,'XMinorTick','on','Box','off');

leg0 = legend('<N_A(t)>_{sim}','<N_B(t)>_{sim}','<N_X(t)>_{sim}');
set(leg0,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');       

%% Plot Extinction time distribution
bins = 10;
figure;
% hist(extt(:,:));                              % histogram for all agents
hist(extt2(:,:),bins);                        % histogram for subspecies

%% Finish
clear temp temp1 temp2 i n x w p1*;
toc
