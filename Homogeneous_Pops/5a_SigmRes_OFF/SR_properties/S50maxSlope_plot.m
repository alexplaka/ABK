%% Make S_50 and gslope vs T (N_TOT) graph

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear;

agents = 0:0.1:200;

k_f = 1;                             % basal forward rate constant
Km_f = 5;                              % MICROSCOPIC MM constant for forward rx
k_r = 0.05;                              % basal reverse rate constant
Km_r = 30;                               % MICROSCOPIC MM constant for reverse rx

% ** Renaming parameters for simplicity **
T = agents;
a = k_f;
b = Km_f;
c = k_r;
d = Km_r;

S_50 = ( c/a .* (T+2*b) ./ (T+2*d) ).^-1;
% dRdS = R*(b+T-R) / ( S*(b+T-2*R) + a/c *(d-T+2*R) );     % General derivative R'(S)
gslope = c/4/a .* T.* (T+2*b).^2 ./ ( b .* (T+2*d) + d .* (T+2*b) );

fig = figure('NumberTitle','off');                          hold on;
set(fig,'Position',[550 1 500 406]);                   
[ax, h1, h2] = plotyy(T,S_50,T,gslope);
set([h1 h2],'LineWidth',2);                       
set(h1,'Color','b','DisplayName','K_{M,f} < K_{M,r}');                          
set(h2,'Color','r');            set(ax(2),'YColor','r');  
xlabel('N_{TOT}');
ylabel(ax(1),'N_{S_{50}}');
ylabel(ax(2),'max dN_R^* / dN_S');
axis(ax(1),[0 T(end) 0 1.1*max(S_50)]);  
set(ax(1),'XMinorTick','on','YMinorTick','on','YTickMode','auto','Box','off');
set(ax(2),'YMinorTick','on','YTickMode','auto','Box','off');

p2 = plot(ax(1),T,a/c*ones(1,size(T,2)),'--','DisplayName','K_{M,f} = K_{M,r}');

leg = legend([h1 p2]);
set(leg,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthWest');


