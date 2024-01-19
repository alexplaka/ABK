function [ON] = detect_SwitchState(all, time, gap, ss, leg_str)
% This function calculates the percentage of simulation runs where
% the (bistable) system is in switch state ON. It then plots the results.
% Input arguments:
% all:      Matrix of size reps x time steps x m, where reps is the number
%           of simulation runs (rows), the number of columns is the number
%           of time steps. The third dimension is m because we are 
%           comparing simulations from m different population structures.
%           Note that all contains the data for a single species.
% time:     Horizontal array of actual time values.
% gap:      It is sometimes useful not to plot all data points in the
%           time series. gap defines how many data points to skip 
%           between plotted points. For fixed time step increments, 
%           this will result in plotting data at dt*gap time unit intervals.
% ss:       Horizontal array with the deterministic fixed points for 
%           switch states ON, OFF. Therefore, for a bistable switch, 
%           the size of ss should 1 x 2.
% leg_str:  String array of size 1 x m, with names of different cases.
%           Used for printing the desired figure legend entries.
%           Ex., leg_str = ["HET1", "HET2"] 
%           Note: string arrays of this format can be used in 
%           Matlab R2017a and later.
% ************************************************************************     
% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

ss_diff = abs( ss(1) - ss(2) );
all_ON = abs( all - ss(1) ) < ss_diff / 4;   % Heuristic cutoff

% Preallocate memory
ON = zeros(size(all,3), size(time,2));
% sdevON = zeros(size(all,3), size(time,2));

fig1 = figure('Name','Mean % ON Trajectories', 'NumberTitle','off');
set(fig1,'Position',[501 1 500 450]);                                   hold on;
set(fig1,'PaperPosition',[0 0 6 5], 'PaperUnits','inches');

for m = 1:size(all, 3)
    
    ON(m,:) = mean(all_ON(:,:,m),1) * 100;
%     sdevON(m,:) = std(all_ON(:,:,m),1) * 100;
    p1(m) = plot(time(1:gap:end), ON(m,1:gap:end),...
        'DisplayName', '$\textrm{' + leg_str(m) + '}$');
%     p2(m) = plot(time(1:gap:end), sdevON(m,1:gap:end), 'DisplayName', '$SDev ' + leg_str(m) + '$');
    
end

xlim([0 time(end)]);                              
set(gca,'XMinorTick','on', 'YMinorTick','on', 'Box','off');
xlabel('$t \; \textrm{(sec)}$', 'Interpreter','LaTeX', 'FontSize',11);                 
ylabel('$\% \; ON$', 'Interpreter','LaTeX','FontSize',11);                             

leg1 = legend(p1); 
set(leg1,'Interpreter','LateX','FontSize',10,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthWest');

% Inner plot
ax2 = axes('Parent',fig1,'Position',[0.62 0.2 0.25 0.25],'FontSize',7);         hold on;
plot(time,ON(1,:) - ON(2,:),'Color','r','Parent',ax2);  
plot([0, time(end)],[0, 0],'--k','LineWidth',1.5);

xlim([0 time(end)]);     ylim([-4, 4]);
xticks([0, 1000, 3000, 5000]);          % adjust accordingly
set(ax2,'XMinorTick','on', 'YMinorTick','on', 'Box','on');
xlabel(ax2,'$t \; \textrm{(sec)}$', 'Interpreter','LaTeX', 'FontSize',8);                 
ylabel(ax2,'$\Delta \% \, ON \; \textrm{(HET1 - HET2)}$',...
    'Interpreter','LaTeX','FontSize',8);                             



