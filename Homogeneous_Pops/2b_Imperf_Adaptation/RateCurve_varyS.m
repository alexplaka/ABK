% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear;

agents = 100;

% Microscopic rates [units: 1/sec]
k1 = 0.50;                       % Synthesis of R; 0th order (upregulated by S; inhibited by X) 
k2 = 0.10;                      % Degradation of R; 1st order (wrt R)
k3 = 0.10;                       % Synthesis of X; 0th order (upregulated by S)
k4 = 0.05;                       % Degradation of X; 1st order 

K = 40;             % K_50 for X inhibiting production of R (# of X to produce half-maximal effect) 
nH = 2;             % Hill coefficient for X inhibiting the production of R

S_level = 0:10:60;

R_array = 0:agents;

% Construct Rate curve: Initialize figure
fig1 = figure('Name','Rate Curve','NumberTitle','off');
set(fig1,'Position',[1 1 500 406]);                         hold on;

for i=1:size(S_level,2)
    
    S = S_level(i);

    % Rate curve must be done for a particular value of X. 
    % This is an arbitrarily chosen value for X, yet one 
    % reflecting its overall dependence on S.
    % Here, I use the steady state value of X for a given value of S.
    X = k3/k4 * S;

    dRdt_degr = k2 * R_array;               
    dRdt_synt = k1 * S * 1 /(1 + (X/K)^nH) * ones(1,size(R_array,2));

    [R_ss_int(i), ssr(i)] = intersections(R_array,dRdt_degr,R_array,dRdt_synt);
    % R_ss_int: Steady-state value of R derived from intersections in Rate curve
    % ssr: Steady-state rate of forward and back reactions

    plot(R_array,dRdt_synt,'Color',[0 0.5 0]);
    plot(R_array,dRdt_degr,'m');
    plot(R_ss_int,ssr,'ok');
    
end

% Finish plot
axis([0 agents 0 35]);         
xlabel('N_R');              ylabel('dN_R / dt');              
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');                   hold off;    
leg1 = legend('Synthesis','Degradation','Location','NorthWest');
set(leg1,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95]);

%% Plot Annotations
% Create textbox
% annotation(fig1,'textbox',...
%     [0.13 0.290 0.12 0.052],'String',{['N_S = ' num2str(S_level(1))]},...
%     'FontWeight','demi','FontName','Times New Roman',...
%     'FitBoxToText','off','LineStyle','none','Color',[0 0.75 0]);

% Create textbox
% annotation(fig1,'textbox',...
%     [0.13 0.407 0.12 0.052],'String',{['N_S = ' num2str(S_level(2))]},...
%     'FontWeight','demi','FontName','Times New Roman',...
%     'FitBoxToText','off','LineStyle','none','Color',[0 0.75 0]);

% Create textbox
% annotation(fig1,'textbox',...
%     [0.13 0.524 0.12 0.052],'String',{['N_S = ' num2str(S_level(3))]},...
%     'FontWeight','demi','FontName','Times New Roman',...
%     'FitBoxToText','off','LineStyle','none','Color',[0 0.75 0]);

% Create textbox
% annotation(fig1,'textbox',...
%     [0.13 0.640 0.12 0.052],'String',{['N_S = ' num2str(S_level(4))]},...
%     'FontWeight','demi','FontName','Times New Roman',...
%     'FitBoxToText','off','LineStyle','none','Color',[0 0.75 0]);


