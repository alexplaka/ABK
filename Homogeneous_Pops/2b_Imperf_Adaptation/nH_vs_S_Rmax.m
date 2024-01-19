% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear;

k1 = 0.50;          % MICROscopic rate of synthesis of R; 0th order [units: 1/sec]
k2 = 0.10;          % MICROscopic rate of degradation of R; 1st order [units: 1/sec]
k3 = 0.10;          % MICROscopic rate of synthesis of X; 0th order [units: 1/sec]
k4 = 0.05;          % Degradation of X; 1st order [units: 1/sec]
K = 40;             % K_50 for X inhibiting production of R (# of X to produce half-maximal effect)

nH_array = 0:0.1:100;

for i=1:size(nH_array,2)
    nH = nH_array(i);
    
    if nH > 1
        S_Rmax(i) = 1 / (nH-1)^(1/nH) / (k3 / k4 / K);
        R_max(i) = (k1 / k2) * S_Rmax(i) / (1 + (k3 / k4 / K)^nH * S_Rmax(i)^nH); 
%         pm(i) = plot(S_Rmax(i),R_max(i),'o','Color',[1 0.75 0.2],...
%             'MarkerFaceColor',[0.87 0.49 0],'DisplayName','max(N_R^*)');
    else
        S_Rmax(i) = NaN;
        R_max(i) = NaN;
    end
 
end

%% Plot nH vs S_Rmax
plot(nH_array,S_Rmax);
axis([0 max(nH_array) 14 max(S_Rmax)]);
set(gca,'XMinorTick','on','YMinorTick','off','Box','off');   
xlabel('n_H');
ylabel('N_{S \rightarrow R^*_{max}}');

% Find nH value where the minimum of S_Rmax occurs
min_ind = find(S_Rmax == min(S_Rmax));
disp(['Minimum occurs at nH = ' num2str(nH_array(min_ind)) ...
    ' , S_Rmax = ' num2str(S_Rmax(min_ind))]);

%% Plot nH vs R_max
plot(nH_array,R_max);
xlabel('n_H');
ylabel('R^*_{max}');
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');   

%% Plot 3D S_Rmax vs R_max vs nH (x vs y vs z)
plot3(S_Rmax,R_max,nH_array,'o','Color',[1 0.75 0.2],'MarkerFaceColor',[0.87 0.49 0]);
xlabel('N_{S \rightarrow R^*_{max}}');
ylabel('N_{R^*_{max}}');
zlabel('n_H');

