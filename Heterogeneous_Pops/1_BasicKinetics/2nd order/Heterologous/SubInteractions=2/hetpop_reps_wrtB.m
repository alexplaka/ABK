%  A + B --> C
% Simulating 2nd order kinetics using my agent-based algorithm.
% The differential rate law for this rx (wrt B being limiting) is used
% This is used to calculate the probability of rx for time duration dt.
% The experiment is repeated in order to obtain average trajectories.

% Heterogeneous Population! Assuming 2 subpopulations.
% Here, I evaluate each possible interaction wrt B.
% (ie, I draw a random number for each possible B-A pair).

% Implementation detail: I randomize the order of sampled B-A interactions.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear; clc;     tic
rng(0);

agents = 10;

maxTime = 200;
dt = 1/100;                                     % Fixed time step increment
t_steps = ceil(maxTime / dt);
    
reps = 1000;                                    % Repeat experiment this # of times

% Initial number of reactant molecules; assume B is limiting.
Ao = agents;                    Bo = ceil(0.7*agents);                
Co = 0;                         Cmax = Bo;

% **************** 2 Distinct Subinteractions *********************
% **********   Set up k-matrix wrt B (size: Bo x Ao)   ************
k_sub = [0.004 , 0.006];                % units: 1/sec
subints = size(k_sub,2);                % number of subinteractions
k1 = k_sub(1)*ones(Bo,Ao/2);            k2 = k_sub(2)*ones(Bo,Ao/2);
k = [k1, k2];
ko = mean(mean(k));                     % units: 1/sec
% *****************************************************************

% ******************** # Subints = Bo*Ao **************************
% **************** Normally-Distributed k *************************
% ********* Set up k-matrix wrt B (size: Bo x Ao) *****************
% ko = 0.001;                     % Mean k value
% sigma = 0.0001;                 % standard deviation of k
% 
% repeat = 1;
% while repeat == 1
%     k = ko + sigma*randn(Bo,Ao);    % Normally distributed with mean=ko and std=sigma
%     repeat = 0;
%     negk = find(k<=0);
%     if isempty(negk)==0             % Ensure there are no negative k values
%         repeat = 1;
%     end
% end
% *****************************************************************

extt1 = zeros(reps, Bo);     % To monitor extinction time of 1st subinteraction
extt2 = zeros(reps, Bo);     % To monitor extinction time of 2nd subinteraction

%% Preallocate memory; Initialize time-dependent sum of molecule numbers
At = zeros(reps,t_steps+1);        Bt = zeros(reps,t_steps+1);        Ct = zeros(reps,t_steps+1);

for n=1:reps

    fprintf(1,'.');
    
    % Arrays for tracking reactant and product agent status
    A = ones(1,Ao);                 B = ones(1,Bo);                 C = zeros(1,Cmax); 
    At(n,1) = sum(A);               Bt(n,1) = sum(B);               Ct(n,1) = sum(C);
    time = zeros(1,t_steps);

    for t=1:t_steps

        tempA = find(A==1);                 % Find 'alive' agents A
        tempB = find(B==1);                 % Find 'alive' agents B
        
        tempA = RandArray(tempA);           % Randomize order of alive A agents
        
        for i = 1:size(tempB,2)
                for j = 1:size(tempA,2)     % Then see if it reacts with an A agent
                        
                        P = 1 - exp(- k(tempB(i),tempA(j)) * dt);       % P_int
%                         P = k(tempB(i),tempA(j)) * dt;                % P_dif
                        
                        if rand < P
                            
                            B(tempB(i)) = 0;       
                            A(tempA(j)) = 0;
                            C(tempB(i)) = 1;   
                            
                            % Assume two subinteractions:
                            if tempA(j) <= Ao/2                 % Extinction time for
                                extt1(n,tempB(i)) = time(t)+dt; % 1st subinteraction
                            else                                % Extinction time for
                                extt2(n,tempB(i)) = time(t)+dt; % 2nd subinteraction
                            end
                            
                            break;    % B(tempB(i)) agent has reacted in this time step.
                            
                        end                                               
                end
        end 

        
    At(n,t+1) = sum(A);         Bt(n,t+1) = sum(B);         Ct(n,t+1) = sum(C);
    time(t+1) = time(t) + dt;
        
    end

end                                         % end 'for reps' loop

% Population statistics of time trajectories
avgA = mean(At);            avgB = mean(Bt);            avgC = mean(Ct);
sdevA = std(At);            sdevB = std(Bt);            sdevC = std(Ct);

% Solve ODE for 2nd order kinetics
[t_sol,y_sol] = ode45(@(t,y) o2_dif(t,y,ko),0:maxTime/100:maxTime,[Ao ; Bo ; Co]);

%% ** Group extinction times by subspecies ** - Works for subspecies = 2 simulation 
% Remove zeros from extt (zeros represent agents that didn't die within totalTime)
f = size(extt1,2);                      
for z=1:f
    extt1_zeros = find(extt1(:,z)==0);
    for a = 1:size(extt1_zeros,1)
        extt1(extt1_zeros(a),z) = NaN;
    end
end

g = size(extt2,2);
for q=1:g
    extt2_zeros = find(extt2(:,q)==0);
    for b = 1:size(extt2_zeros,1)
        extt2(extt2_zeros(b),q) = NaN;
    end
end

extt1 = reshape(extt1,[],1);            % Reshape to column vector
extt2 = reshape(extt2,[],1);
%% Plot time trajectories
fig0 = figure('Name','2nd Order Rx Time course','NumberTitle','off'); 
fig0.PaperUnits = 'inches';
fig0.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.

plot(time,avgA,time,avgB);                                              hold on;
% plot(time,avgC);

% Plot DE Time Courses
% plot(t,E)
plot(t_sol,y_sol(:,1),':c','LineWidth',2);               % scatter(t_sol,y_sol(:,1),'.b');                  
plot(t_sol,y_sol(:,2),':m','LineWidth',2);               % scatter(t_sol,y_sol(:,2),'.g');
% plot(t_sol,y_sol(:,3),':k','LineWidth',2);               % scatter(t_sol,y_sol(:,3),'.r');

p1_devA1 = plot(time,avgA+sdevA,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);
p1_devA0 = plot(time,avgA-sdevA,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);
p1_devB1 = plot(time,avgB+sdevB,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);
p1_devB0 = plot(time,avgB-sdevB,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);
% p1_devC1 = plot(time,avgC+sdevC,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);
% p1_devC0 = plot(time,avgC-sdevC,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);

% title(['N_{A,i} = ' num2str(Ao) ' , N_{B,i} = ' num2str(Bo) ' , N_{X,i} = ' num2str(Co)],...
%     'FontName','Times New Roman','FontSize',12,'FontWeight','Normal');
title(['k_{BA-1} = ' num2str(k_sub(1)) ' sec^{-1} , k_{BA-2} = ' num2str(k_sub(2)) ' sec^{-1}'],...
    'FontName','Times New Roman','FontSize',12,'FontWeight','Normal');

xlabel('t (sec)');                  ylabel('N(t)');          hold off;
axis([0 maxTime 0 agents]);
set(gca,'XMinorTick','on','Box','off');

leg0 = legend('ABK <N_A(t)>','ABK <N_B(t)>','DE N_A(t)','DE N_B(t)');
set(leg0,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');       

% Create textbox
annotation(fig0,'textbox',[0.02 0.94 0.05 0.07],'String',{'b)'},'FontSize',12,...
    'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');

%% Plot Extinction time distribution
bins = 20;
fig1 = figure;
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.

% hist(extt1(:,:),bins);                          % histogram for subinteraction 1
% hist(extt2(:,:),bins);                          % histogram for subinteraction 2

[num1, c1] = hist(extt1,bins);                  % histogram data for subinteraction 1
[num2, c2] = hist(extt2,bins);                  % histogram data for subinteraction 2
Num1 = num1 / (reps*Bo);                        % Normalize
Num2 = num2 / (reps*Bo);                        % Normalize

timeD = 0:maxTime/bins:maxTime;
[t_sol1,y_sol1] = ode45(@(t,y) o2_dif(t,y,k_sub(1)),timeD,[Ao ; Bo ; Co]);
[t_sol2,y_sol2] = ode45(@(t,y) o2_dif(t,y,k_sub(2)),timeD,[Ao ; Bo ; Co]);
% dA = abs(diff(y_sol2(:,1))) / Ao;                     
dB1 = abs(diff(y_sol1(:,2))) / Bo / 2;
dB2 = abs(diff(y_sol2(:,2))) / Bo / 2;

timeDm = maxTime/bins/2:maxTime/bins:maxTime;

b = bar(timeDm,[Num1' Num2']);                                  hold on;
set(b(1),'FaceColor','flat','EdgeColor',[0.70 0.78 1]);         % For subinteraction 1
set(b(2),'FaceColor',[1 0.1 0.1],'EdgeColor',[1 0.6 0.6]);      % For subinteraction 2

plot(timeDm,dB1,'xb');
plot(timeDm,dB2,'or');                                           hold off;

axis([0 maxTime 0 0.3]);
set(gca,'Ytick',0:0.05:0.3,...
    'YtickLabel',{'0' '0.05' '0.10' '0.15' '0.20' '0.25' '0.30'});
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('t (sec) ');             ylabel('Fraction of A+B\rightarrowX Transitions');

leg1 = legend(['Subinteraction BA-1, k_{BA-1} = ' num2str(k_sub(1)) ' sec^{-1}'],...
    ['Subinteraction BA-2, k_{BA-2} = ' num2str(k_sub(2)) ' sec^{-1}']);
set(leg1,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');       

% Coefficient of Determination R^2 for histogram fit to distribution PDF.
disp(['Subint-1: R^2 = ' num2str(CoefDet(Num1,dB1'))]);
disp(['Subint-2: R^2 = ' num2str(CoefDet(Num2,dB2'))]);
%
err = [Num1' - dB1 , Num2' - dB2];              % calculate Deviation

ax2 = axes('Parent',fig1,'Position',[0.62 0.45 0.25 0.25],'FontSize',8);
p_res = plot(timeDm,err,'.','Parent',ax2);                              hold on;
set(p_res(2),'Color','r');
plot(timeDm,zeros(1,size(timeDm,2)),'--','Color',[0.5 0.5 0.5],'Parent',ax2);

axis([0 maxTime -0.04 0.04]);
set(gca,'XMinorTick','on','YMinorTick','on','Box','on');
xlabel('t (sec) ');             ylabel('Deviation');

% Create textbox
annotation(fig1,'textbox',[0.02 0.94 0.05 0.07],'String',{'d)'},'FontSize',12,...
    'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');

%% Finish
clear a b z f temp* g i j n x w p* leg* ax* fig*;
toc
