% This script makes a graph of the scaling factor F as a function 
% of the number of agents X

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear; 

% Parameters
agents = 120;
alpha = [0 0.5 1 2 3];
n = [0.5 1 2];
K = 60;

% figure;             hold on;       % For scatter plot in loop

%% Generate data points
for x = 1:size(alpha,2)
    for y = 1:size(n,2)
        for B=0:agents
            F(x,y,B+1) = (K^n(y) + alpha(x) * B^n(y)) / (K^n(y) + B^n(y));
        end
%         scatter(1:agents,k(x,y,:),'.');
    end
end

%% Plot curves
a0n2  = squeeze(F(1,3,:));
a05n1 = squeeze(F(2,2,:));
a1n   = squeeze(F(3,1,:));
a2n05 = squeeze(F(4,1,:));
a3n2  = squeeze(F(5,3,:));

figure1 = figure('Name','H function','NumberTitle','off');
set(figure1,'Position',[1 1 500 450]);   
set(gcf,'PaperPosition',[0 0 6 5],'PaperUnits','inches');

hold on;

plot(ones(1,10)*K,linspace(0,3,10),'--','Color',[0.8 0.8 0.8]);  % Plot K-50 (vert)
text(K-2,2.75,'$K_{50}$','Interpreter','LaTeX','Rotation',90);

p1 = plot(0:agents,a0n2,'b','LineWidth',2,'DisplayName','\alpha = 0  , \, n_H = 2');
plot(0:K,ones(1,K+1)*a0n2(K+1),'--','Color',[0.8 0.8 0.8]);
a1 = annotation('textbox',[0.735 0.182 0.2 0.047]);
set(a1,'String','$\alpha = 0, \, n_H = 2$','Interpreter','LaTeX',...
    'FontSize',10,'LineStyle','none');

p2 = plot(0:agents,a05n1,'r','LineWidth',2,'DisplayName','\alpha = 0.5, \, n_H = 1');
plot(0:K,ones(1,K+1)*a05n1(K+1),'--','Color',[0.8 0.8 0.8]);
a2 = annotation('textbox',[0.735 0.295 0.22 0.047]);
set(a2,'String','$\alpha = 0.5, \, n_H = 1$','Interpreter','LaTeX',...
    'FontSize',10,'LineStyle','none');

p3 = plot(0:agents,a1n,'g','LineWidth',2,'DisplayName','\alpha = 1');
plot(0:K,ones(1,K+1)*a1n(K+1),'--','Color',[0.8 0.8 0.8]);
a3 = annotation('textbox',[0.80 0.38 0.096 0.047]);
set(a3,'String','$\alpha = 1$','Interpreter','LaTeX',...
    'FontSize',10,'LineStyle','none');

p4 = plot(0:agents,a2n05,'m','LineWidth',2,'DisplayName','\alpha = 2  , \, n_H = 0.5');
plot(0:K,ones(1,K+1)*a2n05(K+1),'--','Color',[0.8 0.8 0.8]);
a4 = annotation('textbox',[0.735 0.54 0.22 0.047]);
set(a4,'String','$\alpha = 2, \, n_H = 0.5$','Interpreter','LaTeX',...
    'FontSize',10,'LineStyle','none');

p5 = plot(0:agents,a3n2,'k','LineWidth',2,'DisplayName','\alpha = 3  , \, n_H = 2');
plot(0:K,ones(1,K+1)*a3n2(K+1),'--','Color',[0.8 0.8 0.8]);
a5 = annotation('textbox',[0.735 0.815 0.2 0.047]);
set(a5,'String','$\alpha = 3, \, n_H = 2$','Interpreter','LaTeX',...
    'FontSize',10,'LineStyle','none');

hold off;

axis([0 agents 0 3]);
set(gca,'YTickLabel',{'0','0.5','1.0','1.5','2.0','2.5','3.0'});
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('$N_X$','Interpreter','LaTeX');                  
ylabel('$F$','Interpreter','LaTeX');

% leg = legend([p1 p2 p3 p4 p5]);
% set(leg,'Location','NorthWest');
% set(leg,'Interpreter','LaTeX','FontSize',9,...
%     'EdgeColor',[0.95 0.95 0.95]);

%% Finish
clear p* a*;
