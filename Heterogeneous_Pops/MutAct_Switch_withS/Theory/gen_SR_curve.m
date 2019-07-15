%   Ep  <--->  E            Rates: k_f, k_r (Reverse rx activated by R)
%    |   ^                  Michaelis-Menten constants: Km_f, Km_r
%    |    \
%    \    |
%    -->  R   -->                      R: Response
%   k_b   ^   k_d
%         | k_s
%         S                            S: Signal

% Simulating mutual activation switch process for R: 
% k_b: R synthesis, 1st order process wrt Ep (but Ep is not consummed)
% k_b: ALTERNATIVE IMPLEMENTATION: R synthesis, 0th order process, UPregulated by Ep
% k_d: R degradation, 1st order process wrt R
% k_s: R synthesis, 1st order process wrt S
% k_s: ALTERNATIVE IMPLEMENTATION: R synthesis, 0th order process, UPregulated by S
% k_f: E synthesis, MM process
% k_r: Ep synthesis, MM process wrt R (but R is not consummed)
% Km_f, Km_r: MM constants for forward and reverse rxs

% Note: Assume number of S molecules/agents is NOT changing (eg, S is DNA molecules)

% Here, we generate the theoretical SR curve for this motif.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

%% Declare variables and functions
clear;                          clc; 

global agents k_b k_d k_s k_f k_r Km_f Km_r S;

agents = 100;
k_b = 0.02;                        % 0th order R synthesis rate, UPregulated by Ep
k_d = 0.075;                       % 1st order degradation rate (units: 1/sec)
k_s = 0.05;                        % 0th order R synthesis rate UPregulated by S
k_f = 1;                           % basal forward rate (1st order)
k_r = 0.05;                        % basal reverse rate (1st order)
Km_f = 5;                          % MICROSCOPIC Michaelis-Menten constant for forward rx
Km_r = 10;                         % MICROSCOPIC Michaelis-Menten constant for reverse rx

S_array = 0:0.1:30;

R_ss = zeros(size(S_array,2),3);
Ep_ss = zeros(size(S_array,2),3);

for s=1:size(S_array,2)
    
    S = S_array(s);
    
    [R_ss(s,:), Ep_ss(s,:)] = RateCurve_s;
        
end             % end 'for s' loop
%% Remove duplicate FPs from monostable systems
bistab = [];
for i=1:size(R_ss,1)
   if R_ss(i,1) == R_ss(i,2) && isempty(bistab)==1
       R_ss(i,2:3) = NaN;
       Ep_ss(i,2:3) = NaN;
   elseif R_ss(i,1) == R_ss(i,2) && isempty(bistab)==0
       R_ss(i,1:2) = NaN;
       Ep_ss(i,1:2) = NaN;
   else
       bistab = [bistab i];
   end
end

disp(['Bistability at S = ' num2str(S_array(bistab(1))) ' : ' num2str(S_array(bistab(end)))]);

%% Plot
fig1 = figure('Name','SR curve','NumberTitle','off');
set(fig1,'Position',[501 1 500 406]);                       hold on;

for h=1:size(R_ss,2)
    if mod(h,2) == 0
        p(h) = plot(S_array,R_ss(:,h),'--','LineWidth',2,...
            'Color',[1 0.60 0.39],'DisplayName','Unstable');
    else
    p(h) = plot(S_array,R_ss(:,h),'-b','LineWidth',2,'DisplayName','Stable');
    end
end

axis([S_array(1) S_array(end) 0 max(max(R_ss))]);        hold off;
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('N_S');         ylabel('N_R^*');         

leg1 = legend([p(1) p(2)]);
set(leg1,'Location','SouthEast');
set(leg1,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95]);

%% Finish
clear fig* leg* p* i h s