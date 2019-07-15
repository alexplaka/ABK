% Compare different forms of computing the transition probability.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

clear;clc;

dt = 0:0.002:0.1;
k = 0.01;

figure('Position',[0 0 1100 400]);
%% First run
At = 1000;

for i=1:size(dt,2)
    P_can(i) = 1 - 1 / (1 + k * (At-1) * dt(i));
    P_int(i) = 1 - 1 /( At * (1-exp(-k*dt(i))) + exp(-k*dt(i)) );
    P_ber(i) = 1 - exp(-k * (At-1) * dt(i));
    P_dif(i) = k * (At-1) * dt(i);
end

subplot(1,2,1);             
plot(dt,P_can,'m',dt,P_int,'b',dt,P_ber,'g',dt,P_dif,'r');                          hold on;
plot(0.01*ones(1,21),0.02:0.05:1.02,'LineStyle','--','Color',[0.8 0.8 0.8]);        hold off;
axis([0 dt(end) 0 1]);              
title(['k_2 = ' num2str(k) ' sec^{-1} ; N_A(t_n) = ' num2str(At)],'FontName','Times New Roman','FontSize',12);
xlabel('\Deltat (sec)');                  ylabel('P_{2A\rightarrowX}');
% leg = legend('P_{can}','P_{int}','P_{ber}','P_{dif}','Location','NorthWest');
% set(leg,'FontName','Times New Roman');        set(leg,'FontSize',9);
set(gca,'Box','off');   
set(gca,'XTick',[0:0.01:0.10]);     set(gca,'YTick',[0:0.1:1]);
set(gca,'XTickLabel',{'0','0.01','0.02','0.03','0.04','0.05','0.06','0.07','0.08','0.09','0.10'});      
set(gca,'YTickLabel',{'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'});      
set(gca,'XMinorTick','on');         set(gca,'YMinorTick','on');

%% Second run
At = 10000;

for i=1:size(dt,2)
    P_can(i) = 1 - 1 / (1 + k * (At-1) * dt(i));
    P_int(i) = 1 - 1 /( At * (1-exp(-k*dt(i))) + exp(-k*dt(i)) );
    P_ber(i) = 1 - exp(-k * (At-1) * dt(i));
    P_dif(i) = k * (At-1) * dt(i);
end

subplot(1,2,2);                     
plot(dt,P_can,'m',dt,P_int,'b',dt,P_ber,'g',dt,P_dif,'r');                          hold on;
plot(0.01*ones(1,21),0.02:0.05:1.02,'LineStyle','--','Color',[0.8 0.8 0.8]);        hold off;
axis([0 dt(end) 0 1]);              
title(['k_2 = ' num2str(k) ' sec^{-1} ; N_A(t_n) = ' num2str(At)],'FontName','Times New Roman','FontSize',12);
xlabel('\Deltat (sec)');                  ylabel('P_{2A\rightarrowX}');
leg = legend('P_{can}','P_{int}','P_{ber}','P_{dif}','Location','SouthEast');
set(leg,'FontName','Times New Roman');        set(leg,'FontSize',9);
set(gca,'Box','off');
set(gca,'XTick',[0:0.01:0.10]);     set(gca,'YTick',[0:0.1:1]);
set(gca,'XTickLabel',{'0','0.01','0.02','0.03','0.04','0.05','0.06','0.07','0.08','0.09','0.10'});      
set(gca,'YTickLabel',{'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'});      
set(gca,'XMinorTick','on');         set(gca,'YMinorTick','on');
