%   Ep  <--->  E            Rates: k_f, k_r (Reverse rx activated by R)
%    |   ^                  
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
% k_f: E synthesis
% k_r: Ep synthesis (but R is not consummed)

% Note: Assume number of S molecules/agents is NOT changing (eg, S is DNA molecules)

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

%% Declare variables and functions
clear;                          clc; 

global agents k_b k_d k_s k_f k_r S;

agents = 100;
k_b = 0.01;                        % 0th order R synthesis rate, UPregulated by Ep
k_d = 0.01;                       % 1st order degradation rate (units: 1/sec)
k_s = 0.01;                        % 0th order R synthesis rate UPregulated by S
k_f = 0.1;                           % basal forward rate (1st order)
k_r = 0.05;                        % basal reverse rate (1st order)

S_array = 0:0.1:100;

R_ss = zeros(1,size(S_array,2));

for s=1:size(S_array,2)
    
    S = S_array(s);
    
    % Clumped parameters for determining R_ss
    a = 1;
    b = - k_s/k_d * S + k_f/k_r - k_b/k_d * agents;
    c = - k_f*k_s / (k_r*k_d);

    R_ss_exact(1) = (-b + sqrt(b^2 - 4*a*c)) / (2*a);
    R_ss_exact(2) = (-b - sqrt(b^2 - 4*a*c)) / (2*a);
    R_ss(s) = R_ss_exact(find(R_ss_exact>0));
        
end             % end 'for s' loop
%% Plot
fig1 = figure('Name','SR curve','NumberTitle','off');
set(fig1,'Position',[501 1 500 406]);                       hold on;

plot(S_array,R_ss,'-b','LineWidth',2);

axis([S_array(1) S_array(end) 0 max(max(R_ss))]);        hold off;
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('N_S');         ylabel('N_R^*');         
%% Finish
clear fig* leg* p* i h s