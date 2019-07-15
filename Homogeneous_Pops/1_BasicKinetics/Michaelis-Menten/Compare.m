% Here, I compare the ODE solutions and ABK time trajectories of the 
% Simple and Full implementations of Michaelis-Menten kinetics.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;

global k_f k_r k_cat;
global Km E_tot;
agents = 20;
E_tot = 5;

Sa(1) = agents;                     Sb(1) = 0;
Se(1) = E_tot;                      Sea(1) = 0;

k_f = 0.01;                         k_r = 0.05;
k_cat = 0.05;                       
Km = (k_r + k_cat) / k_f;

%% Compare DE solutions
finaltime = 300;
cd('Simple Treatment');
[t_sim, y_sim] = ode45(@mm_sim_dif,0:finaltime/100:finaltime,[Sa(1) ; Sb(1)]);
cd('../Full Treatment');
[t_full, y_full] = ode45(@mm_full_dif,0:finaltime/100:finaltime,[Sa(1) ; Sb(1); Se(1); Sea(1)]);
cd ..;

%% Plot DE solutions
figure1 = figure('Name','Solution to DE(s)','NumberTitle','off');
set(figure1,'Position',[1 1 500 450]);                                              hold on;
plot(t_full,y_full(:,1),'-b');              plot(t_full,y_full(:,2),'-r'); 
plot(t_sim,y_sim(:,1),'--c');                plot(t_sim,y_sim(:,2),'--','Color',[1 0.6 0]);
xlabel('t (sec)');                          ylabel('N(t)');                         hold off;
axis([0 finaltime 0 Sa(1)]);                % axis tight;    
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
% title('Solution to Differential Equation(s)','FontName','Times New Roman','FontSize',11)
leg = legend('N_A(t) full','N_B(t) full','N_A(t) abbr', 'N_B(t) abbr');      
set(leg,'Location','East');
set(leg,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95]);

%% Compare average ABK solutions of time trajectories: Get data first
cd('Full Treatment');                       
load(['Ao=' num2str(Sa(1)) '_reps.mat'],'avgA','avgB','time');   % Make sure file exists
avgA_full = avgA;                           avgB_full = avgB;
time_full = time;
cd('../Simple Treatment');                  
load(['Ao=' num2str(Sa(1)) '_reps.mat'],'avgA','avgB','time');   % Make sure file exists
avgA_sim = avgA;                            avgB_sim = avgB;
time_sim = time;
cd ..;                                      clear avgA avgB time;

%% Plot ABK time trajectories
figure2 = figure('Name','ABK Results','NumberTitle','off');
set(figure2,'Position',[1 1 500 450]);                                              hold on;
plot(time_full,avgA_full,'-b');             plot(time_full,avgB_full,'-r');
plot(time_sim,avgA_sim,'--c');              plot(time_sim,avgB_sim,'--','Color',[1 0.6 0]);
xlabel('t (sec)');                          ylabel('N(t)');                         hold off;
axis([0 finaltime 0 Sa(1)]);                % axis tight;             
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
% title('ABK algorithm','FontName','Times New Roman','FontSize',11)
leg = legend('<N_A(t)>_{sim} Full','<N_B(t)>_{sim} Full',...
    '<N_A(t)>_{sim} Abbr', '<N_B(t)>_{sim} Abbr');      
set(leg,'Location','East');
set(leg,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95]);
%% Plot both DE and ABK trajectories
figure;                                                                             hold on;
% Plot ABK trajectories first
plot(time_sim,avgA_sim,'-k');               plot(time_sim,avgB_sim,'-g');
plot(time_full,avgA_full,'-b');             plot(time_full,avgB_full,'-r');

% Plot DE trajectories
plot(t_sim,y_sim(:,1),':k');                plot(t_sim,y_sim(:,2),':g');
plot(t_full,y_full(:,1),':b');              plot(t_full,y_full(:,2),':r');

xlabel('t (sec)');                          ylabel('N(t)');                         hold off;
axis([0 finaltime 0 Sa(1)]);                % axis tight;             
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
leg = legend('N_A(t) simple','N_B(t) simple','N_A(t) full', 'N_B(t) full');      
set(leg,'Location','East');
set(leg,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95]);
