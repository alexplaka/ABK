%% Plot [selected] time course with detected peaks shown
% Using Report Nomenclature: X = R , Rp = Xp  (script name = report name) in legend.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

trial = 20;                         % Choose simulation run/trial (trial <= reps)

fig0 = figure('Name','Time Course','NumberTitle','off');                    hold on;
set(fig0,'Position',[1 1 500 450]);
set(fig0,'PaperPosition',[0 0 6 5],'PaperUnits','inches');

plot(time,X_all(trial,:));
plot(time,Rp_all(trial,:));
% findpeaks(X_all(trial,:),time,'MinPeakProminence',30,'MinPeakHeight',70,'MinPeakDistance',45);
% findpeaks(Rp_all(trial,:),time,'MinPeakProminence',8,'MinPeakHeight',15,'MinPeakDistance',45);

% *** Smooth out trajectories ***
smoothTrajX = sgolayfilt(X_all(trial,:),5,6001);
smoothTrajRp = sgolayfilt(Rp_all(trial,:),5,6001);
% plot(time,smoothTrajX,'c','Linewidth',2);
% plot(time,smoothTrajRp,'k','Linewidth',2);
findpeaks(smoothTrajX,time,'MinPeakProminence',26,'MinPeakHeight',65,'MinPeakDistance',45);
findpeaks(smoothTrajRp,time,'MinPeakProminence',7,'MinPeakHeight',14,'MinPeakDistance',45);
% *******************************

% axis([0 time(end) 0 N_max]);  
axis([0 time(end) 0 165]);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off','XGrid','off','YGrid','off');
xlabel('$t \; \textrm{(sec)}$','Interpreter','LaTeX');                 
ylabel('$N(t)$','Interpreter','LaTeX');                                       hold off;

title(['$\textbf{Normally-distributed delay: } \bar{\delta} =' num2str(lag_mean) ...
    '\, \textrm{sec} \, \textrm{,} \; \sigma_{\delta} =' num2str(lag_std) ...
    '\, \textrm{sec} $'],'Interpreter','LaTeX');

leg0 = legend('$\textrm{ABK } N_R(t)$','$\textrm{ABK } N_{Xp}(t)$',...
    '$\textrm{Smoothed ABK } N_R(t)$','$N_R \textrm{ detected peaks}$',...
    '$\textrm{Smoothed ABK } N_{Xp}(t)$','$N_{Xp} \textrm{ detected peaks}$',...
    'Location','NorthEast');
% leg0 = legend('ABK Sample Run N_R(t)','ABK Sample Run N_{Xp}(t)','Location','Best'); % Without DEs
set(leg0,'Interpreter','LaTeX','FontSize',10,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');

%% Plot histogram of peak intervals (of specified 'trial')
fig1 = figure('Name','Time Course','NumberTitle','off');                    hold on;
set(fig1,'Position',[1 1 500 450],'PaperPosition',[0 0 6 5],'PaperUnits','inches');

% [pksX, locsX] = findpeaks(X_all(trial,:),time,...
%     'MinPeakProminence',40,'MinPeakHeight',80,'MinPeakDistance',50);
% [pksRp, locsRp] = findpeaks(Rp_all(trial,:),time,...
%     'MinPeakProminence',10,'MinPeakHeight',20,'MinPeakDistance',50);

[pksX, locsX] = findpeaks(smoothTrajX,time,...
    'MinPeakProminence',40,'MinPeakHeight',75,'MinPeakDistance',50);
[pksRp, locsRp] = findpeaks(smoothTrajRp,time,...
    'MinPeakProminence',10,'MinPeakHeight',20,'MinPeakDistance',50);


peakIntervalX = diff(locsX);
peakIntervalRp = diff(locsRp);

[numX, cX] = hist(peakIntervalX);
[numRp, cRp] = hist(peakIntervalRp);

bar(cX,numX);       
bar(cRp,numRp,'r');