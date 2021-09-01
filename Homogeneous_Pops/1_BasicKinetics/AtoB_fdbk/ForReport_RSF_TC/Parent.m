% Make a plot consisting of several data sets (different alpha values)

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;   clc;           tic
rng(0);
pb = waitbar(0,'0');

global alpha;
a = [0 1 2];

for x = 1:size(a,2)
    
    progress = x/size(a,2);
    waitbar(progress,pb,sprintf('%.0f%%',progress*100));

    alpha = a(x);
    reps_child;
    
    trajTime = time;
    trajA(x,:) = avgA;
    trajAsdev(x,:) = sdevA;
    DE_t = t_sol;
    DE_A(x,:) = y_sol(:,1);
end

%% Plot Time Courses
colors = ['rbgmkc'];
figure1 = figure('Name','o1 rx w/ Feedback ','NumberTitle','off');
set(figure1,'Position',[1 1 500 450]);                  hold on;

for x = 1:size(a,2)
    p(x) = plot(trajTime,trajA(x,:),colors(x),'MarkerSize',3,'DisplayName',['\alpha = ' num2str(a(x))]);                 hold on;
%     p_dev1 = plot(trajTime,trajA(x,:)+trajAsdev(x,:),'LineStyle','--','Color',[0.8 0.8 0.8]);
%     p_dev0 = plot(trajTime,trajA(x,:)-trajAsdev(x,:),'LineStyle','--','Color',[0.8 0.8 0.8]);

    pd = plot(DE_t,DE_A(x,:),':','Color',colors(x));        

end

xlabel('t (sec)');                 ylabel('N_A(t)');                 hold off;
title(['n_H = ' num2str(n)],'FontName','Times New Roman','FontSize',11);
axis tight;                         % axis([0 t_max 0 agents]); 
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
leg = legend([p(:)]);
set(leg,'Location','NorthEast');
set(leg,'FontName','Times New Roman','FontSize',10,...
    'EdgeColor',[0.95 0.95 0.95]);

%% Finish
close(pb);
toc