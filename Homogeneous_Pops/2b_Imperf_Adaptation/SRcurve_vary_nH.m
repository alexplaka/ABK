%    k3     k4
%    -->  X -->                        
%   /    /
%   S   /                               S: Signal
%   \  /
%    -->  R -->                         R: Response
%    k1     k2

% Simulating birth-death process for R and X, and X inhibits the production of R: 
% k1: synthesis of R, 0th order process (upregulated by S; inhibited by X)
% k2: degradation of R, 1st order process (R)
% k3: synthesis of X, 0th order process (upregulated by S)
% k4: degradation of X, 1st order process (X)

% Note: Assume number of S molecules/agents is NOT changing
% Note: Assume number of X molecules/agents is NOT changing while inhibiting formation of R

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;     % clc;    

color = 'brgkcmy'; 

k1 = 0.50;          % MICROscopic rate of synthesis of R; 0th order [units: 1/sec]
k2 = 0.10;          % MICROscopic rate of degradation of R; 1st order [units: 1/sec]
k3 = 0.10;          % MICROscopic rate of synthesis of X; 0th order [units: 1/sec]
k4 = 0.05;          % Degradation of X; 1st order [units: 1/sec]
K = 40;             % K_50 for X inhibiting production of R (# of X to produce half-maximal effect)

S = 0:50;
nH_array = 0:1:4;

R_ss = zeros(size(nH_array,2),size(S,2));
S_Rmax = zeros(1,size(nH_array,2));
R_max = zeros(1,size(nH_array,2));

fig1 = figure('Name','SR Curve','NumberTitle','off');
set(fig1,'Position',[1 1 500 406]);                                         hold on;

for i=1:size(nH_array,2)
    
    nH = nH_array(i);             % Hill coefficient for X inhibiting the production of R
    R_ss(i,:) = (k1 / k2) .* S ./ (1 + (k3 / k4 / K)^nH .* S.^nH); 
    
    p(i) = plot(S,R_ss(i,:),'Color',color(i),'LineWidth',2,...
        'DisplayName',['n_H = ' num2str(nH)]);
    
    
end

% Redo of above loop to make sure maximum of nH=1 curve is "above" other curves (Matlab bug?)
for i=1:size(nH_array,2)
    nH = nH_array(i);
    
    if nH > 1
        S_Rmax(i) = 1 / (nH-1)^(1/nH) / (k3 / k4 / K);
        R_max(i) = (k1 / k2) * S_Rmax(i) / (1 + (k3 / k4 / K)^nH * S_Rmax(i)^nH); 
        pm(i) = plot(S_Rmax(i),R_max(i),'o','Color',[1 0.75 0.2],...
            'MarkerFaceColor',[0.87 0.49 0],...
            'DisplayName','max N_R^*');
    else
        S_Rmax(i) = NaN;
        R_max(i) = NaN;
    end
 
end

%% Draw how max N_R* changes as nH changes
nH_array2 = 1.1:0.4:100;

for i=1:size(nH_array2,2)
    nH = nH_array2(i);
    S_Rmax2(i) = 1 / (nH-1)^(1/nH) / (k3 / k4 / K);
    R_max2(i) = (k1 / k2) * S_Rmax2(i) / (1 + (k3 / k4 / K)^nH * S_Rmax2(i)^nH); 
end
pm2 = plot(S_Rmax2,R_max2,'--','Color',[1 0.75 0.2],...
    'MarkerFaceColor',[0.87 0.49 0]);

%% Finish plot
axis([0 S(end) 0 100]);
xlabel('N_S');              ylabel('N_R^*');              
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');                   hold off;    
leg1 = legend([p pm(end)],'Location','NorthWest');
set(leg1,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95],...
    'Interpreter','tex');

clear i j p* fig* leg*;