%  k_b   k_d
%  --> R -->                        R: Response
%      ^
%      | k_s
%      S                            S: Signal

% Simulating birth-death process for R, and its synthesis from S: 
% k_b: birth, 0th order process; k_d: death, 1st order process
% k_s: synthesis from S, 1st order process
% Note: Assume number of S molecules/agents is NOT changing (eg, S is DNA molecules)
% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;      clc;     tic;
rng(0);
pb = waitbar(0,'0');

agents = 1200;
reps = 100;

N_avo = 6.02e23;                % Avogadro's number
V = 1e-21;                      % Volume in L

totalTime = 700;                % Simulation time (sec)
dt = 1/50;                      % Fixed time step increment (sec)

if exist('dt','var') == 1         
    t_steps = totalTime / dt;
else
    t_steps = 1200;             % empirically chosen value (for lambda = 0.8)
end
time = zeros(1,t_steps);

global k_b k_d k_s S;
k_bm = 3.322e-5;                % Molar "birth" rate; 0th order (units: M/sec)
k_b = k_bm * N_avo * V;         % Microscopic 0th order birth rate
k_d = 0.02;                     % "Death" rate; 1st order (units: 1/sec)
k_s = 1;                        % Synthesis rate; 1st order (units: 1/sec)

S_array = 0:20;                 % Do a simulation for each of these values of S
a = size(S_array,2);

% Preallocate memory
R_ss = zeros(1,a);              ABK_r = zeros(reps,a);

for n = 1:a
    
    progress = n / a;
    waitbar(progress,pb,sprintf('%.0f%% %',progress*100));
    
    S = S_array(n);             % Assume number of S molecules/agents is NOT changing

    % Deterministic Steady-State value for R
    R_ss(n) = (k_b + k_s * S) / k_d;

    %   ---------- Time-Independent Probability values -------------  
    P_b = k_b * dt;                 % P_dif of "birth" process (0th order)

    P_d = 1 -  exp(-k_d * dt);      % P_ber of "death" process (1st order wrt R)
    % P_d = k_d * dt;                 % P_dif of "death" process (1st order wrt R)

    P_s0 = k_s * S * dt;            % P_dif of synthesis process (0th order wrt S)
    %   ------------------------------------------------------------  

    for w=1:reps
        
        Tr = zeros(1,t_steps);              
        Tr(1) = 0;                  % ** Initial condition for R **

    %   ---------------------- Set up SRV --------------------------  
        tempR = zeros(1,agents);             
        % Put "1" where agents are "alive", then randomize the array
        for c=1:Tr(1),                  tempR(c)=1;             end
        tempR = RandArray(tempR);       % Randomize array
        % Markov process, so only previous and current time steps needed --> 2 rows:    
        R = [tempR ; tempR];            % Initialize vector storing state of A agents
        clear c tempR;
    %   ------------------------------------------------------------  

        for t = 2:t_steps

            if rand < P_b                           % "Birth", 0th order reaction
                temp = find(R(1,:)==0);             % Randomly choose R agent which becomes "alive" 
                R(2,temp(ceil(rand * size(temp,2)))) = 1;
            end

            if rand < P_s0                          % "Birth", 0th order reaction promoted by S
                temp = find(R(1,:)==0);             % Randomly choose R agent which becomes "alive" 
                R(2,temp(ceil(rand * size(temp,2)))) = 1;
            end    

            for i = 1:size(R,2)
                if R(1,i) == 1                      % if R agent still exists          
                    if rand < P_d                   % "Death", 1st order reaction
                        R(2,i) = 0;                 % R agent is degraded
                    end
                end
            end

            Tr(t) = sum(R(2,:));
            R(1,:) = R(2,:);
            time(t) = time(t-1) + dt;
        end

    %   Find average R of last 300sec (steady-state has been reached)
        ABK_r(w,n) = mean(Tr(end-300*50:end));

    end             % end "for w" loop
end             % end "for n" loop

avg_R = mean(ABK_r);                        std_R = std(ABK_r);   
%% Determine fit of obtained SR curve to theoretical expectation 
[m, CC] = FitCurve_lineb_s(S_array,avg_R,k_b/k_d);   % m = slope, CC = CorrCoeff
percError =  abs(m - k_s/k_d) / (k_s/k_d);
disp(['Corr Coeff R = ' num2str(CC)]);
disp(['% error(m) = ' num2str(percError)]);

%% Plot SR curve
figure1 = figure('Name','Signal-Response Curve','NumberTitle','off');
set(figure1,'Position',[1 1 500 406]);                           hold on;

plot(S_array,R_ss,S_array,avg_R,'.r');
plot(S_array,avg_R+std_R,'--','Color',[0.8 0.8 0.8]);
plot(S_array,avg_R-std_R,'--','Color',[0.8 0.8 0.8]);

xlabel('N_S');    ylabel('N_R^*'); 

axis([0 S_array(end) 0 1000]);
% axis([0 10 0 500]);

set(gca,'XTick',0:2:20);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
leg = legend('DE','ABK <N_R^*>');
set(leg,'Location','NorthWest');
set(leg,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95]);

%% Finish
close(pb);
clear progress pb a h i n t temp S;
toc
