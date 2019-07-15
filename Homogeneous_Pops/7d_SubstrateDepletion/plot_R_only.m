% Run this script after '*_reps.m' has been run.

%% Plot [selected] time course for R only
trial = 4;

figure('Name',['Time Course, S=' num2str(S)],'NumberTitle','off','Position',[1 1 500 406]);     
hold on;

p_Rs = plot(time(1:20:end),R_all(trial,1:20:end),'r','DisplayName','ABK N_R(t)');                                

p_Rd = plot(t_sol,y_sol(:,1),'--','Color',[0.75 0 1],'LineWidth',2,'DisplayName','DE N_R(t)');                                       

axis tight;                                     % axis([0 finaltime 0 agents]);                               
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');                           hold off;
xlabel('t (sec)');                 ylabel('N(t)'); 

leg0 = legend([p_Rs , p_Rd]);       % Just R
set(leg0,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthWest');                                              

%% Plot AVG trajectory for R
figure('Name',['Time Course, S=' num2str(S)],'NumberTitle','off','Position',[1 1 500 406]);           
hold on;

plot(time(1:20:end),avg_R(1:20:end),'r','LineWidth',2);
plot(t_sol,y_sol(:,1),'--','Color',[0.75 0 1],'LineWidth',2);   
plot(time(1:20:end),avg_R(1:20:end)+sdev_R(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
% plot(time(1:20:end),avg_R(1:20:end)-sdev_R(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);

% axis tight;  
axis([0 time(end) 0 170]);

set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('t (sec)');            ylabel('N_R(t)');                                        

leg1 = legend('<N_R(t)>_{sim}','DE  N_R(t)','Location','NorthEast');
set(leg1,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95]);                                              hold off;
