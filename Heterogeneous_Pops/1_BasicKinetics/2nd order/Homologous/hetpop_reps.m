%  A + A --> C
% Simulating 2nd order kinetics using my agent-based algorithm.
% The differential rate law for this rx is used
% This is used to calculate the probability of rx for time duration dt.
% The experiment is repeated in order to obtain average trajectories

% Heterogeneous Population!
% Here, I evaluate each possible interaction wrt A.
% (ie, I draw a random number for each possible A-A pair).

% Extinction times are logged and a histogram is constructed.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear; clc;     tic
rng(0);

global DE_canon;
agents = 8;

maxTime = 200;
dt = 1/100;                            % Fixed time step increment
t_steps = ceil(maxTime / dt);

reps = 2000;                              % Repeat experiment this # of times

% Initial number of reactant molecules; assume B is limiting.
Ao = agents;                              Co = 0;
Cmax = Ao;

% **************** 2 Distinct Subinteractions ****************
% ********* Set up k-matrix wrt A (size: Ao x Ao) ************
k_sub = [0.004 , 0.006];
[k , abundance] = gen_randPHM(Ao,k_sub);

subints = size(k_sub,2); 
extt1 = zeros(reps, Ao);     % To monitor extinction time of A agents
extt2 = zeros(reps, Ao);     % To monitor extinction time of A agents
% ************************************************************

% **************** # Subinteractions = Agents ****************
% **************** Normally-Distributed k ********************
% ********* Set up k-matrix wrt A (size: Ao x Ao) ************
% ko = 0.005;                     % Mean k value
% sigma = 0.000;                 % standard deviation of k
% 
% repeat = 1;
% while repeat == 1
%     k = ko + sigma*randn(Ao);       % Normally distributed with mean=ko and std=sigma
%     repeat = 0;
%     negk = find(k<=0);
%     if isempty(negk)==0             % Ensure there are no negative k values
%         repeat = 1;
%     end
% end
% 
% subints = Ao * (Ao - 1); 
% extt = zeros(reps,Ao);              % To monitor extinction time of A agents
% ************************************************************
ko = mean(mean(k));          % Bimolecular microscopic kinetic constant [units: 1/sec];

% Preallocate memory; Initialize time-dependent sum of molecule numbers
At = zeros(reps,t_steps+1);                Ct = zeros(reps,t_steps+1);

%% Run simulation

for n=1:reps

    fprintf(1,'.');
    
    % Arrays for tracking reactant and product agent status
    A = ones(1,Ao);                 C = zeros(1,Cmax);  % Assume Co = 0;
    At(n,1) = sum(A);               Ct(n,1) = sum(C);
    time = zeros(1,t_steps);

    for t=1:t_steps

        for i = 1:size(A,2)
            if A(i) == 1                    % if this A agent is alive
                for j = 1:size(A,2)           % then see if it reacts with another A agent
                    if i ~= j && A(j) == 1            % if this other A agent is alive
                        
                        P = 1 - exp(- k(i,j) * dt);
                        
                        if rand < P
                            
                            A(i) = 0;       
                            A(j) = 0;
                            C(j) = 1;  
                            
                            if subints == Ao*(Ao-1)      % Normally-distributed k      
                                extt(n,i) = time(t)+dt;      % Extinction time
                                extt(n,j) = time(t)+dt;      % Extinction time
                                
                            elseif subints == 2     % For only 2 distinct subinteractions 
                                if k(i,j) == k_sub(1);            % Extinction time for
                                    extt1(n,i) = time(t)+dt;      % 1st subinteraction
                                else                              % Extinction time for
                                    extt2(n,i) = time(t)+dt;      % 2nd subinteraction
                                end
                            end
                            
                            
                            break;   % Make sure A(i) is not reevaluated in this time step
                            
                        end
                    end
                end
            end
        end 
        
    At(n,t+1) = sum(A);                  Ct(n,t+1) = sum(C);
    time(t+1) = time(t) + dt;
        
    end

end                 % end "for reps" loop

% Population statistics of time trajectories
avgA = mean(At);                        avgC = mean(Ct);
sdevA = std(At);                        sdevC = std(Ct);

%% ** Group extinction times by subspecies ** - Works for subinteractions = 2 simulation 
if subints == Ao*(Ao-1)
    for z=1:agents
        extt_zeros = find(extt(:,z)==0);
        for a = 1:size(extt_zeros,1)
            extt(extt_zeros(a),z) = NaN;
        end
    end
elseif subints == 2
% Remove zeros from extt1/2 (zeros represent agents that didn't die within totalTime)
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
end

%% Solve ODE for 2nd order kinetics
% DE_canon = 0;               % Agent-based DE for the process
% [t_sol0,y_sol0] = ode45(@(t,y) o2_dif(t,y,ko),0:maxTime/100:maxTime,[Ao ; Co]);

DE_canon = 1;               % Canonical DE for the process
[t_sol1,y_sol1] = ode45(@(t,y) o2_dif(t,y,ko),0:maxTime/100:maxTime,[Ao ; Co]);

%% Plot time trajectories
fig0 = figure('Name','2nd Order Rx Time course','NumberTitle','off'); 
fig0.PaperUnits = 'inches';
fig0.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.

plot(time,avgA);                                              hold on;
% plot(time,avgC);

% Plot DE Time Courses
% plot(t_sol0,y_sol0(:,1),':k','LineWidth',2);                                
plot(t_sol1,y_sol1(:,1),':c','LineWidth',2);               

p1_devA1 = plot(time,avgA+sdevA,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);
p1_devA0 = plot(time,avgA-sdevA,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);
% p1_devC1 = plot(time,avgC+sdevC,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);
% p1_devC0 = plot(time,avgC-sdevC,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);

% title(['N_{A,i} = ' num2str(Ao) ' , N_{X,i} = ' num2str(Co)],...
%     'FontName','Times New Roman','FontSize',12,'FontWeight','Normal');
% title(['k_{BA-1} = ' num2str(k_sub(1)) ' sec^{-1} , k_{BA-2} = ' num2str(k_sub(2)) ' sec^{-1}'],...
%     'FontName','Times New Roman','FontSize',12,'FontWeight','Normal');

xlabel('t (sec)');                  ylabel('N(t)');          hold off;
axis([0 maxTime 0 agents]);
set(gca,'XMinorTick','on','Box','off');

leg0 = legend('<N_A(t)>_{sim}','DE N_A(t)');                                  % Just species A
% leg0 = legend('<N_A(t)>_{sim}','<N_X(t)>_{sim}','DE N_A(t)','DE N_X(t)');   % All species
set(leg0,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');       

% Create textbox
annotation(fig0,'textbox',[0.02 0.94 0.05 0.07],'String',{'a)'},'FontSize',12,...
    'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');

%% Plot Extinction time distribution
bins = 20;
fig1 = figure;
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.

if subints == Ao*(Ao-1)
    hist(extt(:,:),bins);                         % histogram for all agents
    
elseif subints == 2
    % hist(extt1(:,:),bins);                          % histogram for subinteraction 1
    % hist(extt2(:,:),bins);                          % histogram for subinteraction 2

    [num1, c1] = hist(extt1,bins);                  % histogram data for subinteraction 1
    [num2, c2] = hist(extt2,bins);                  % histogram data for subinteraction 2
    Num1 = num1 / (reps*Ao);                        % Normalize
    Num2 = num2 / (reps*Ao);                        % Normalize
% ** Problem with Normalization! sum(Num1) or (Num2) do not approach 1 as maxTime increases. **
% *********************

    timeD = 0:maxTime/bins:maxTime;
%     DE_canon = 0;
    [t_sol1p,y_sol1p] = ode45(@(t,y) o2_dif(t,y,k_sub(1)),timeD,[Ao ; Co]);
    [t_sol2p,y_sol2p] = ode45(@(t,y) o2_dif(t,y,k_sub(2)),timeD,[Ao ; Co]);
    dA1 = abs(diff(y_sol1p(:,2))) / Ao * abundance(1);      % Normalize
    dA2 = abs(diff(y_sol2p(:,2))) / Ao * abundance(2);      % Normalize

    timeDm = maxTime/bins/2:maxTime/bins:maxTime;

    b = bar(timeDm,[Num1' Num2']);                                  hold on;
    set(b(1),'FaceColor','flat','EdgeColor',[0.70 0.78 1]);         % For subinteraction 1
    set(b(2),'FaceColor',[1 0.1 0.1],'EdgeColor',[1 0.6 0.6]);      % For subinteraction 2

    plot(timeDm,dA1,'b');
    plot(timeDm,dA2,'r');                                           hold off;

    axis([0 maxTime 0 0.20]);
    set(gca,'Ytick',0:0.05:0.3,...
        'YtickLabel',{'0' '0.05' '0.10' '0.15' '0.20' '0.25' '0.30'});
    set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
    xlabel('t (sec) ');             ylabel('Fraction of 2A\rightarrowX Transitions');

    leg1 = legend(['Subinteraction AA-1, k_{AA-1} = ' num2str(k_sub(1)) ' sec^{-1}'],...
        ['Subinteraction AA-2, k_{AA-2} = ' num2str(k_sub(2)) ' sec^{-1}']);
    set(leg1,'FontName','Times New Roman','FontSize',9,...
        'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');       

    % Coefficient of Determination R^2 for histogram fit to distribution PDF.
    disp(['Subint-1: R^2 = ' num2str(CoefDet(Num1,dA1'))]);
    disp(['Subint-2: R^2 = ' num2str(CoefDet(Num2,dA2'))]);
    
    err = [Num1' - dA1 , Num2' - dA2];              % calculate Deviation

    ax2 = axes('Parent',fig1,'Position',[0.62 0.45 0.25 0.25],'FontSize',8);
    p_res = plot(timeDm,err,'.','Parent',ax2);                              hold on;
    set(p_res(2),'Color','r');
    plot(timeDm,zeros(1,size(timeDm,2)),'--','Color',[0.5 0.5 0.5],'Parent',ax2);

    axis([0 maxTime -0.02 0.02]);
    set(gca,'XMinorTick','on','YMinorTick','on','Box','on');
    xlabel('t (sec) ');             ylabel('Deviation');

    % Create textbox
    annotation(fig1,'textbox',[0.02 0.94 0.05 0.07],'String',{'d)'},'FontSize',12,...
        'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');
end

%% Finish
clear a b z f temp* g i j n x w p* leg* ax* fig* P *zeros *_sol* A DE_canon;
toc

% Result: It works? Yes!
% This implementation gives the same average trajectory of A as the 
% the simulation of the homogeneous population.
% For heterogeneous populations, obtained similar results to A+B->X case.
% Nore problem in normalization of PDF in histogram.