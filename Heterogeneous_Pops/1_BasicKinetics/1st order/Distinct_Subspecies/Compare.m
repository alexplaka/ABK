% Compare SDev and CV results

clear; clc;   

dirs = 5;
AVG = zeros(dirs,1001);                     AVG_TH = zeros(dirs,1001);  
SDEV = zeros(dirs,1001);                    SDEV_TH = zeros(dirs,1001);
CV = zeros(dirs,1001);                      CV = zeros(dirs,1001);
CV_poisson = zeros(dirs,1001);

k_val = ['0.1-0.9' ; '0.2-0.8' ; '0.3-0.7' ;'0.4-0.6' ; '0.5-0.5'];

for a=1:dirs
    cd(['reps=2000_k=' k_val(a,:)]);
    load('A=10_subsp=2.mat','k','time','avg','sdev','cv','agents','totalTime');
    cd ..;
    [y_theory_avg , y_theory_sdev] = Het_o1_Predict(time,agents,k); 
    AVG(a,:) = avg;                 AVG_TH(a,:) = y_theory_avg;
    SDEV(a,:) = sdev;               SDEV_TH(a,:) = y_theory_sdev;
    CV(a,:) = cv;
    CV_poisson(a,:) = 1./sqrt(avg);
    CV_TH = SDEV_TH ./ AVG_TH;
end
clear avg sdev cv a;
%% Compare simulation standard deviation 
figure3 = figure('Name','SDev Comparison','NumberTitle','off');
set(figure3,'Position',[501 501 500 450]);                        hold on;    
figure3.PaperUnits = 'inches';
figure3.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.

for b=1:dirs-1
    p3(b) = plot(time,SDEV(b,:)/agents,'LineWidth',2,'DisplayName',...
        ['$k_{A1} = ' k_val(b,1:3) '\, \textrm{sec}^{-1} ; k_{A2} = ' k_val(b,5:7) '\, \textrm{sec}^{-1} $']);                          
    p3_th(b) = plot(time,SDEV_TH(b,:)/agents,'--','LineWidth',1,...
        'DisplayName','$\textrm{Prediction}$','Color',[0.8 0.2 0.2]);
end

% Last plot is for a homogeneous population
p3(dirs) = plot(time,SDEV(dirs,:)/agents,'LineWidth',2,'Color','c','DisplayName',...
        ['$\bar{k} = ' k_val(dirs,1:3) '\, \textrm{sec}^{-1} \; \textrm{(Homogeneous pop.)}$']);                          
p3_th(dirs) = plot(time,SDEV_TH(dirs,:)/agents,'--','LineWidth',1,...
        'DisplayName','$\textrm{Prediction}$','Color',[0.8 0.2 0.2]);

xlabel('$t \; \textrm{(sec)}$','Interpreter','LaTeX');           
ylabel('$SDev \Big( < \!N_A(t) \! > \Big) \, / \, N_{A,i}$','Interpreter','LaTeX');                 
set(gca,'XTick',0:totalTime,'XMinorTick','on','YMinorTick','on','Box','off');
set(gca,'YTick',0:0.02:0.18);
set(gca,'YTickLabel',{'0','0.02','0.04','0.06','0.08','0.10','0.12','0.14','0.16','0.18'});
axis([0 totalTime 0 0.18]);                                        hold off;

leg3 = legend([p3(end:-1:1) p3_th(dirs)]);
set(leg3,'Interpreter','LaTeX','EdgeColor',[0.95 0.95 0.95],...
    'Location','NorthEast','FontSize',10)

% Create textbox
annotation(figure3,'textbox',...
    [0.00771 0.93471 0.05010 0.06429],'String','a)',...
    'FontSize',12,'FontName','Times New Roman',...
    'FitBoxToText','off','EdgeColor','none');

%% Plot coefficient of variation
figure4 = figure('Name','CV Comparison','NumberTitle','off');
set(figure4,'Position',[501 501 500 450]);                        hold on;    
figure4.PaperUnits = 'inches';
figure4.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.

% colors = ['c','m','k','g'];

for c=1:dirs-1
    p4(c) = plot(time,CV(c,:),'-.','LineWidth',2,...
        'DisplayName',...
        ['$k_{A1} = ' k_val(c,1:3) '\, \textrm{sec}^{-1} ; k_{A2} = ' k_val(c,5:7) '\, \textrm{sec}^{-1} $']);               
%     pt4_th(c) = plot(time,CV_TH(c,:),':',...
%         'Color',[0.6 0.6 0.6],'LineWidth',2);    
%     pt4_pois(c) = plot(time,CV_poisson(c,:),'--','Color',[0.2 0.75 0.5],...
%         'DisplayName','$\textrm{Poisson}$','LineWidth',1);
end

% Last plot is for a homogeneous population
p4(dirs) = plot(time,CV(dirs,:),'-.','LineWidth',2,'Color','c',...
        'DisplayName',['$\bar{k} = ' k_val(dirs,1:3) '\, \textrm{sec}^{-1} \; \textrm{(Homogeneous pop.)}$']);               
% pt4_th(dirs) = plot(time,CV_TH(dirs,:),':','Color',[0.6 0.6 0.6],...
%         'DisplayName','$\textrm{Prediction}$','LineWidth',2);
% pt4_pois(dirs) = plot(time,CV_poisson(dirs,:),'--','Color',[0.2 0.75 0.5],...
%         'DisplayName','$\textrm{Poisson}$','LineWidth',1);
      
axis([0 totalTime 0 4]);        
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
set(gca,'XTick',0:1:totalTime);
set(gca,'YTick',0:0.5:4,'YTickLabel',{'0','0.5','1.0','1.5','2.0','2.5','3.0','3.5','4.0'});

leg4 = legend([p4(end:-1:1)]); %#ok<NBRAK>
set(leg4,'Interpreter','LaTeX','FontSize',10,'EdgeColor',[0.95 0.95 0.95]);
set(leg4,'Location','NorthWest');

xlabel('$t \; \textrm{(sec)}$','Interpreter','LaTeX');           
ylabel('$\eta$','Interpreter','LaTeX','FontSize',11);                                                   hold off;                                    

% Create textbox
annotation(figure4,'textbox',...
    [0.00771 0.93471 0.05010 0.06429],'String','b)',...
    'FontSize',12,'FontName','Times New Roman',...
    'FitBoxToText','off','EdgeColor','none');
