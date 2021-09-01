% H + S --> S + S                       Rate constant k_c  (2nd order)
%     S --> D                           Rate constant k_d  (1st order)
% Simulating Kermack-McKendrick model using the agent-based algorithm (ABK).

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear; clc;  tic
rng(0);
pb = waitbar(0,'0');

global k_c k_d;
k_c = 0.03;                 % k_c: 2nd order rate constant [units: 1/sec], c = Contagion
k_d = 0.05;                 % k_d: 1st order rate constant [units: 1/sec], d = Death

t_max = 50;                             % in seconds
dt = 1/100;                             % Fixed time step increment
steps = t_max / dt;

reps = 500;                               % Repeat simulation this number of times
% Preallocate species-dependent time trajectory Sums/matrices (rows = single run time trajectories)
Sh = zeros(reps,steps);     Ss = zeros(reps,steps);     Sd = zeros(reps,steps);    

% Initial number of Healthy, Sick, Dead people
Ho = 15;                   So = 1;                     Do = 0;
agents = Ho + So + Do;                  % Assuming 'agents' stays constant

b = k_d / (k_c * Ho);                   % See Strogatz Problem 3.7.6 (2nd ed)
if b < 1
    disp('Epidemic predicted.');
else
    disp('No Epidemic predicted.');
end

for n=1:reps

    progress = n/reps;
    waitbar(progress,pb,sprintf('%.0f%% %',progress*100));
        
    % Initialize time-dependent sum of agents
    Ht = zeros(1,steps);        St = zeros(1,steps);        Dt = zeros(1,steps);
    Ht(1) = Ho;                 St(1) = So;                 Dt(1) = Do;

    tempH = zeros(1,agents);    tempS = zeros(1,agents);    tempD = zeros(1,agents);              
    % Put "1" where agents are "alive", then randomize the arrays. Make sure arrays are complementary!
    for h=1:Ho,              tempH(h)=1;         end;            tempH = RandArray(tempH);
    notH = RandArray(find(tempH == 0));
    for s=1:So,              tempS(notH(s))=1;   end;            % tempS = RandArray(tempS);
    notHS = find(tempH==0 & tempS==0);
    for d=1:Do,              tempD(notHS(d))=1;  end;            % tempD = RandArray(tempD);                     
    Hv = [tempH ; tempH];               % Initialize vector storing state of H agents  
    Sv = [tempS ; tempS];               % Initialize vector storing state of S agents
    Dv = [tempD ; tempD];               % Initialize vector storing state of D agents
    % - Markov process, so only previous and current time steps needed --> 2 rows:          
    clear h s d tempH tempS tempD notH notHS;

    time = zeros(1,steps);

%   Time-Independent Probability values with respect to S
%     P_d = k_d * dt;                             % P_dif-S of 1st order S --> D
    P_d = 1 - exp(-k_d * dt);                   % P_ber-S of 1st order S --> D
    
    for t=2:steps                   % while Ht(t-1) > 0
                
%       Time-dependent Probability values with respect to S
        P_c(t-1) = k_c * Ht(t-1) * dt;               % P_dif-S of 2nd order H + S --> S + S
%         P_c(t-1) = 1 - exp(-k_c * Ht(t-1) * dt);     % P_ber-S of 2nd order H + S --> S + S
           
        tempH = find(Hv(1,:)==1);      tempS = find(Sv(1,:)==1);    % tempD = find(Dv(1,:)==1);
        
        for i = 1:size(tempS,2)
            r = rand;
            if r < P_c(t-1)                             % 2nd order H + S --> S + S   
                h = ceil(rand*size(tempH,2));           % Pick random Healthy agent
                if h==0,        continue;         end
                Hv(2,tempH(h)) = 0;                     % Healthy agent ...
                Sv(2,tempH(h)) = 1;                     % becomes Sick
                Sv(2,tempS(i)) = 1;                     % Sick agent stays alive, but sick
            elseif r >= P_c(t-1) && r < P_c(t-1)+P_d    % 1st order S --> D
                Sv(2,tempS(i)) = 0;                     % Sick agent ...                
                Dv(2,tempS(i)) = 1;                     % becomes Dead
            end 
        end
        
        Ht(t) = sum(Hv(2,:));        St(t) = sum(Sv(2,:));         Dt(t) = sum(Dv(2,:));  
        
        time(t) = time(t-1) + dt;
        Hv(1,:) = Hv(2,:);           Sv(1,:) = Sv(2,:);            Dv(1,:) = Dv(2,:);  
    end
    
    Sh(n,:) = Ht;                Ss(n,:) = St;               Sd(n,:) = Dt; 

end         % end 'for n' loop

avgH = mean(Sh);                     sdevH = std(Sh);
avgS = mean(Ss);                     sdevS = std(Ss);
avgD = mean(Sd);                     sdevD = std(Sd);

%% Solve ODE for mixed order kinetics
[t_sol,y_sol] = ode45(@KM_dif,0:t_max/1000:t_max,[Ho; So; Do]);
%% Plot Average Time Courses (n = reps)
figure1 = figure('Name','Kermack-McKendrick Avg Trajectories','NumberTitle','off');
set(figure1,'Position',[1 1 500 450]);
p1 = plot(time,avgH,'b','MarkerSize',3,'DisplayName','<N_H(t)>_{sim}');                 hold on;
% p1_dev1 = plot(time,avgH+sdevH,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p1_dev0 = plot(time,avgH-sdevH,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p1s = plot(time,Sh(end,:),'m','MarkerSize',3,'DisplayName','Sample Trajectory');

p2 = plot(time,avgS,'r','MarkerSize',3,'DisplayName','<N_S(t)>_{sim}');                 
% p2_dev1 = plot(time,avgS+sdevS,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p2_dev0 = plot(time,avgS-sdevS,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p2s = plot(time,Ss(end,:),'m','MarkerSize',3,'DisplayName','Sample Trajectory');

p3 = plot(time,avgD,'g','MarkerSize',3,'DisplayName','<N_D(t)>_{sim}');                 
% p3_dev1 = plot(time,avgD+sdevD,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p3_dev0 = plot(time,avgD-sdevD,'LineStyle','--','Color',[0.8 0.8 0.8]);
% p3s = plot(time,Sd(end,:),'m','MarkerSize',3,'DisplayName','Sample Trajectory');

p1d = plot(t_sol,y_sol(:,1),'-.b','DisplayName','DE');        
p2d = plot(t_sol,y_sol(:,2),'-.r');
p3d = plot(t_sol,y_sol(:,3),'-.g');

xlabel('t (sec)');                 ylabel('N(t)');                  hold off;   
axis([0 t_max 0 agents]);                          % axis tight;                                        
title(['N_{H,i} = ' num2str(Ho) ' , N_{S,i} = ' num2str(So) ' , N_{D,i} = ' num2str(Do)],...
    'FontName','Times New Roman','FontSize',11)
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
leg = legend([p1 p2 p3]);
% set(leg,'OuterPosition',[0.673 0.374 0.182 0.123]);
set(leg,'Location','NorthWest');
set(leg,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95]);
% set(p??,'Visible','off');
%% Plot all trajectories of a given species
% figure2 = figure('Name','Kermack-McKendrick Trajectories','NumberTitle','off');
% set(figure1,'Position',[501 1 500 450]);                                        hold on;
% 
% for n=1:reps
%     plot(time,Sh(n,:),'b');
% %     plot(time,Ss(n,:),'r');
% %     plot(time,Sd(n,:),'c');
% end
%% Find trajectories that deviate from DE predictions
fewDeaths = find(Sd(:,end) <= 1);
deviation = size(fewDeaths,1) / reps * 100;
%% HS State-space (2D)
figure('Name','State Space','NumberTitle','off');               hold on;
scatter(y_sol(:,1),y_sol(:,2),'.k');       % plot deterministic trajectory
% scatter(avgH,avgS,3,'.c');                 % plot average (n=reps) stochastic trajectory
trial = 1;
scatter(Sh(trial,:),Ss(trial,:),'or');     % plot specific stochastic trajectory
scatter(Sh(trial+5,:),Ss(trial+5,:),'.b'); % plot specific stochastic trajectory
% comet(Sh(trial,:),Ss(trial,:),0.9);
axis([0 agents 0 12]);
% legend('');
xlabel('N_H');                    ylabel('N_S');
%% Finish
close(pb);
clear progress temp* x deltaC deltaD r n h i p* *v Ht St Dt;
toc
