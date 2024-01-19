%  A --> B  ; Simulating 1st order kinetics. 
% Considering individual agents.
% *** Heterogeneous population of A ***

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear; clc;             tic
rng(0);

% Assume constant time step increment (sec)
dt = 0.01;
totalTime = 10;         % Simulation time (sec)
t_steps = totalTime / dt;

reps = 2000;              % Number of times to repeat experiment
% k_calc = zeros(size(dt_array,2),n);
S = zeros(reps,t_steps);

agents = 10;            % Population size
%  # of agents should be integrally divisible by # of subspecies

extt = zeros(reps, agents);     % To monitor extinction time of agents

% ****** Set population's k values ************
% Specify subspecies k values (horizontal array)
k_sub = [0.3 , 0.7];
subspecies = size(k_sub,2);                     % number of subspecies

k = [];
for h=1:subspecies
    k = [k , k_sub(h)*ones(1,agents/subspecies)];
end
ko = mean(k);
% *********************************************

addpath(pwd);
dir_name = ['reps=' num2str(reps) '_k=' num2str(k_sub(1)) '-' num2str(k_sub(2))];
mkdir(dir_name);
cd(dir_name);

% Time-independent probability for a given agent (vectorized)
P = 1 - exp(-k .* dt);                  % P from integrated rate law
% P = k .* dt;                            % P from differential rate law

for n = 1:reps

    A = ones(2,agents);     % Only 2 rows because of Markov property      
    S(n,1) = agents;
    time = zeros(1,t_steps);
    
    for t = 1:t_steps
        
        tempAa = find(A(1,:)==1);               % find "alive" A agents
        
        for i = 1:size(tempAa,2)
            if rand < P(tempAa(i))              % ** check probability condition **
                A(2,tempAa(i)) = 0;             % agent dies
                extt(n,tempAa(i)) = time(t)+dt;                    
            end             % 'if rand' statement
        end                 % 'for i' loop
        
        S(n,t+1) = sum(A(2,:));
        
        A(1,:) = A(2,:);
        time(t+1) = time(t) + dt;
    end                     % 'for t' loop
end                         % 'for n' loop

avg = mean(S);                          
sdev = std(S);
cv = sdev ./ avg;                               % Coefficient of variation

clear A i h n t temp* dir_name;
save('A=10_subsp=2.mat');                      
%% ** Run code from this point forward when loading saved data file! **
%% ** Group extinction times by subspecies **
% Remove zeros from extt (zeros represent agents that didn't die within totalTime)
for z=1:agents
    extt_zeros = find(extt(:,z)==0);
    for a = 1:size(extt_zeros,1)
        extt(extt_zeros(a),z) = NaN;
    end
end

extt1=reshape(extt,[],1);                           % Extinction times for ENTIRE species A

extt2 = zeros(reps*agents/subspecies,subspecies);   % Extinction times by subspecies of A
s = 1;              w = 1;
for j=1:agents
    extt2((w-1)*reps+1:w*reps,s) = extt(:,j);
    
    if j < agents && k(j) ~= k(j+1)
        s = s + 1;
        w = 1;
    else
        w = w + 1;
    end  
end

%% Calculate curves for average k (ko) and subspecies k values
y_theory_ko = agents*exp(- ko * time);
% y_theory_sub1 = agents * exp(- k_sub(1) * time);
% y_theory_sub2 = agents * exp(- k_sub(2) * time);

%% Calculate expected time trajectory (and SDEV) using my theoretical approach:
[y_theory_avg , y_theory_sdev] = Het_o1_Predict(time,agents,k); 
disp(['ABK-Theory time trajectory: R^2 = ' num2str(CoefDet(avg,y_theory_avg))]);
%% Plot average time trajectory
figure1 = figure('Name','Time trajectory','NumberTitle','off');
set(figure1,'Position',[1 1 500 450]);                              hold on;    
figure1.PaperUnits = 'inches';
figure1.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.

p1 = plot(time,avg,'b','LineWidth',2,'DisplayName','$\textrm{ABK } <\!N_A(t)\!>$');                              
plot(time,avg+sdev,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);
plot(time,avg-sdev,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);

p2 = plot(time,y_theory_ko,'g','LineStyle',':','LineWidth',2,...
    'DisplayName','$N_A(t) = 10 \, e^{- \bar{k}\, t}$');
% p3 = plot(time,y_theory_sub1,'-r','MarkerSize',5);
% p4 = plot(time,y_theory_sub2,'-k','MarkerSize',5);

p5 = plot(time(1:10:end),y_theory_avg(1:10:end),':r','LineWidth',2,...
    'DisplayName','Prediction');
% plot(time,y_theory_avg+y_theory_sdev,'r','LineStyle','-','LineWidth',2);
% plot(time,y_theory_avg-y_theory_sdev,'r','LineStyle','-','LineWidth',2);

set(gca,'XTick',0:1:totalTime,'XMinorTick','on','Box','off');
title(['$N_{A,i} = ' num2str(agents) '$'],'Interpreter','LaTeX','FontSize',12);
axis([0 10 0 agents]);
xlabel('$t \; \textrm{(sec)}$','Interpreter','LaTeX');                      
ylabel('$N_A$','Interpreter','LaTeX');                              hold off;

leg1 = legend([p1 p2 p5]);
set(leg1,'Interpreter','LaTeX','FontSize',10,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');       

% Create textbox
annotation(figure1,'textbox',...
    [0.00771 0.93471 0.05010 0.06429],'String','a)',...
    'FontSize',12,'FontName','Times New Roman',...
    'FitBoxToText','off','EdgeColor','none');

hgsave('TC.fig');

% **** Plot sample time trajectories ****
% trial = [1 2];                  
% plot(time,S(trial(1),:),'ob','MarkerSize',3);                       hold on;
% plot(time,S(trial(2),:),'+c','MarkerSize',3);
% xlabel('t (sec)');                      ylabel('N_A');

%% Plot Extinction time distribution for each subspecies in the same plot
bins = 10;

figure2 = figure('Name','Extinction Event Distribution','NumberTitle','off');
set(figure2,'Position',[501 1 500 450]);                            hold on;    
figure2.PaperUnits = 'inches';
figure2.PaperPosition = [0 0 6 5];              % Control the size of printed (eps) fig.

% hist(extt(:,:));                              % histogram for all agents
% hist(extt2(:,:),bins);                        % histogram for subspecies
% h = findobj(gca,'Type','patch');              % use set(h(1),...) commands below to customize

exp_PDF1 = k_sub(1) / (bins/totalTime) * exp(- k_sub(1) * time);
exp_PDF2 = k_sub(2) / (bins/totalTime) * exp(- k_sub(2) * time);
% PDFs need to be adjusted for the number of bins in the histogram
% in order to obtain proper fits. Adjustment: divide by (bins/totalTime)

[num1, c1] = hist(extt2(:,1),bins);             % histogram data for subspecies 1
[num2, c2] = hist(extt2(:,2),bins);             % histogram data for subspecies 2
Num1 = num1 / (reps*agents/subspecies);         % Normalize
Num2 = num2 / (reps*agents/subspecies);         % Normalize

% b = bar([c1' c2'],[Num1' Num2']);         % Works in Matlab 2013a
b = bar(totalTime/bins/2:totalTime/bins:totalTime,[Num1' Num2']);           
set(b(1),'FaceColor',[0.3 0.75 0.93],'EdgeColor',[0.70 0.78 1],...  % For subspecies 1
    'DisplayName',['$\textrm{Subspecies } A1 \textrm{, } k_{A1} = ' num2str(k_sub(1)) '\, \textrm{sec}^{-1}$']);    
set(b(2),'FaceColor',[0.93 0.69 0.13],'EdgeColor',[1 0.6 0.6],...   % For subspecies 2
    'DisplayName',['$\textrm{Subspecies } A2 \textrm{, } k_{A2} = ' num2str(k_sub(2)) '\, \textrm{sec}^{-1}$']);
pdf(1) = plot(time,exp_PDF1,'-','Color',[0.4 0.4 0.4],'LineWidth',1,...
    'DisplayName','$\textrm{Subspecies } A1 \textrm{, PDF}$');
pdf(2) = plot(time,exp_PDF2,'-','Color',[0.84 0.3 0.8],'LineWidth',1,...
    'DisplayName','$\textrm{Subspecies } A2  \textrm{, PDF}$');

% *** Use curve-fitting to verify Exponential Dist PDF. ***
% [est1,y_fit1] = FitCurve_exp_Fixk(c1',Num1',k_sub(1));    % fit for fixed k
% [est2,y_fit2] = FitCurve_exp_Fixk(c2',Num2',k_sub(2));    % fit for fixed k
% [est1,y_fit1] = FitCurve_exp(c1',Num1');
% [est2,y_fit2] = FitCurve_exp(c2',Num2');
% * Graph Fits *
% plot(c1,y_fit1,'-','Color',[0.15 0.15 1],'LineWidth',1);
% plot(c2,y_fit2,'-','Color',[1 0.15 0.15],'LineWidth',1);
% *********************************************************

set(gca,'XTick',[0 1 2 3 4 5 6 7 8 9 10],...
    'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('$t \; \textrm{(sec)}$','Interpreter','LaTeX');                      
ylabel('$\textrm{Fraction of } A \rightarrow X \textrm{ Transitions}$','Interpreter','LaTeX');

leg2 = legend([b pdf]);
set(leg2,'Interpreter','LaTeX','FontSize',10,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');      

% Create textbox
annotation(figure2,'textbox',...
    [0.00771 0.93471 0.05010 0.06429],'String','b)',...
    'FontSize',12,'FontName','Times New Roman',...
    'FitBoxToText','off','EdgeColor','none');

hgsave('Event_Dist.fig');

% Note: A histogram of the extinction events for the entire species
% is unnecessary given the agreement with the expected subspecies
% distributions. However, I have plotted it below. The predicted
% distribution is numerically determined from calling Het_o1_Predict
% (it can also be determined analytically).

% Coefficients of Determination R^2 for histogram fits to exponential dist PDFs.
exp_PDF1d = k_sub(1) / (bins/totalTime) * exp(- k_sub(1) * c1); % Discretize exp PDF
exp_PDF2d = k_sub(2) / (bins/totalTime) * exp(- k_sub(2) * c2); % Discretize exp PDF

disp(['A1 Extinction time distribution: R^2 = ' num2str(CoefDet(Num1,exp_PDF1d))]);
disp(['A2 Extinction time distribution: R^2 = ' num2str(CoefDet(Num2,exp_PDF2d))]);

%% Plot Extinction time distribution for the ENTIRE species A
% bins = 10;
% 
% figure;                                              hold on;
% 
% [A_avg1 , A_sdev1] = Het_o1_Predict(0:0.5:totalTime,agents,k);      %
% dA = abs(diff(A_avg1)) / agents; 
% % PDFs need to be adjusted for the number of bins in the histogram
% % in order to obtain proper fits. Adjustment: divide by (bins/totalTime)
% 
% [num, c] = hist(extt1,bins);                  % histogram data for ENTIRE species A
% Num = num / (reps*agents);                    % Normalize
% 
% b = bar(totalTime/bins/2:totalTime/bins:totalTime,Num');               
% set(b(1),'FaceColor','flat','EdgeColor',[0.70 0.78 1]);         % For species A
% plot(totalTime/bins/2:totalTime/bins:totalTime,dA,'-','Color',[0.55 0.55 1],'LineWidth',2);
% 
% set(gca,'XTick',[0 1 2 3 4 5 6 7 8 9 10],...
%     'XMinorTick','on','YMinorTick','on','Box','off');
% xlabel('t (sec) ');             ylabel('Fraction of A\rightarrowX Transitions');
% 
% % leg1 = legend(['Subspecies A1, k_{A1} = ' num2str(k_sub(1)) ' sec^{-1}'],...
% %     ['Subspecies A2, k_{A2} = ' num2str(k_sub(2)) ' sec^{-1}']);
% % set(leg1,'FontName','Times New Roman','FontSize',9,...
% %     'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');       
% clear A_avg1 A_sdev1;

%% Compare simulation standard deviation with my theoretical prediction
figure3 = figure('Name','SDev Comparison','NumberTitle','off');
set(figure3,'Position',[501 501 500 450]);                        hold on;    
figure3.PaperUnits = 'inches';
figure3.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.

p31(1) = plot(time,sdev/agents,'LineWidth',2,'DisplayName','$\textrm{ABK}$');                          
p31(2) = plot(time,y_theory_sdev/agents,'--r','LineWidth',2,...
    'DisplayName','$\textrm{Prediction}$');

xlabel('$t \; \textrm{(sec)}$','Interpreter','LaTeX');           
ylabel('$SDev \Big( < \!N_A(t) \! > \Big) \, / \, N_{A,i}$','Interpreter','LaTeX');                 
set(gca,'XTick',0:totalTime,'XMinorTick','on','YMinorTick','on','Box','off');
set(gca,'YTick',0:0.02:0.18,'YTickLabel',{'0','0.02','0.04','0.06','0.08','0.10','0.12','0.14','0.16','0.18'});
axis([0 totalTime 0 0.18]);                                        hold off;

leg3 = legend(p31);
set(leg3,'Interpreter','LaTeX','EdgeColor',[0.95 0.95 0.95],...
    'Location','NorthEast','FontSize',10)

% ax2_inset = axes('Parent',figure3,'Position',[0.4 0.19 0.25 0.25],'FontSize',7.5);        hold on;
% p2_inset = plot(time,cv,'k','DisplayName','ABK','LineWidth',2,'Parent',ax2_inset);               
% pt2_inset(1) = plot(time,1./sqrt(avg),':','Color',[0.2 0.75 0.5],...
%     'DisplayName','Poisson','LineWidth',2,'Parent',ax2_inset);
% 
% axis([0 totalTime 0 2.0]);        
% set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
% set(gca,'XTick',0:2:totalTime,'FontSize',7);
% set(gca,'YTick',0:0.4:2,'YTickLabel',{'0','0.4','0.8','1.2','1.6','2.0'});
% leg2_inset = legend([p2_inset pt2_inset(1)]);
% set(leg2_inset,'FontName','Times New Roman','FontSize',9,...
%     'EdgeColor',[0.95 0.95 0.95]);
% set(leg2_inset,'Position',[0.53 0.21 0.18 0.06]);
% % set(leg2_inset,'Location','NorthEast');
% 
% xlabel('$t \; \textrm{(sec)}$','Interpreter','LaTeX');           
% ylabel('$\eta$','Interpreter','LaTeX','FontSize',10);                                                   hold off;                                    

% Create textbox
annotation(figure3,'textbox',...
    [0.00771 0.93471 0.05010 0.06429],'String','c)',...
    'FontSize',12,'FontName','Times New Roman',...
    'FitBoxToText','off','EdgeColor','none');

hgsave('SDev.fig');

% *** Find Coefficient of Determination, R^2 ***
disp(['SDev: R^2 = ' num2str(CoefDet(sdev,y_theory_sdev))]);
% **********************************************

%% Plot coefficient of variation (separate plot)
figure4 = figure('Name','SDev Comparison','NumberTitle','off');
set(figure4,'Position',[501 501 500 450]);                        hold on;    
figure4.PaperUnits = 'inches';
figure4.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.

p4 = plot(time,cv,'k','DisplayName','$\textrm{ABK}$','LineWidth',2);               
pt4(1) = plot(time,1./sqrt(avg),':','Color',[0.2 0.75 0.5],...
    'DisplayName','$\textrm{Poisson}$','LineWidth',2);
% pt4(2) = plot(time,1./sqrt(y_theory_avg),'-','Color',[0.2 0.75 0.5],...
%     'DisplayName','$\textrm{Poisson}$','LineWidth',2);

axis tight;
% axis([0 totalTime 0 2]);        
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
set(gca,'XTick',0:1:totalTime);
% set(gca,'YTick',0:0.4:2,'YTickLabel',{'0','0.4','0.8','1.2','1.6','2.0'});
leg4 = legend([p4 pt4]);
set(leg4,'Interpreter','LaTeX','FontSize',10,'EdgeColor',[0.95 0.95 0.95]);
set(leg4,'Location','SouthEast');

xlabel('$t \; \textrm{(sec)}$','Interpreter','LaTeX');           
ylabel('$\eta$','Interpreter','LaTeX','FontSize',10);                                                   hold off;                                    

% Create textbox
annotation(figure4,'textbox',...
    [0.00771 0.93471 0.05010 0.06429],'String','d)',...
    'FontSize',12,'FontName','Times New Roman',...
    'FitBoxToText','off','EdgeColor','none');

hgsave('CV.fig');

% *** Find Coefficient of Determination, R^2 ***
disp(['CV: R^2 = ' num2str(CoefDet(cv,1./sqrt(avg)))]);
% **********************************************

%% Finish
clear a b h i j s w z leg* p* temp* A;
cd ..;
toc