%  A + B --> C
% Simulating 2nd order kinetics using my agent-based algorithm.
% The differential rate law for this rx (wrt B being limiting) is used
% This is used to calculate the probability of rx for time duration dt.
% The experiment is repeated in order to obtain average trajectories.
% Here, I evaluate each possible interaction wrt B.
% (ie, I draw a random number for each possible B-A pair).

% **** Heterogeneous Population represented as a graph. ****
% **** PHM is a sparse matrix. ****

% Implementation detail: I randomize the order of sampled B-A interactions.  

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear; clc;     tic
rng(0);

agents = 10;

maxTime = 200;
dt = 1/100;                            % Fixed time step increment
t_steps = ceil(maxTime / dt);

reps = 100;                              % Repeat experiment this # of times

% Initial number of reactant molecules; assume B is limiting.
Ao = agents;                    Bo = ceil(0.7*agents);                
Co = 0;                         Cmax = Bo;

% **************** # Subspecies = Agents *********************
% **************** Normally-Distributed k ********************
% ********* Set up k-matrix wrt B (size: Bo x Ao) ************
ko = 0.001;                     % Mean k value
sigma = 0.0001;                 % standard deviation of k
u = [1 2 4 5 5 7 7];            v = [1 9 2 1 5 5 8];        

repeat = 1;
while repeat == 1
    z = ko + sigma*randn(1,size(i,2));     % Normally distributed with mean=ko and std=sigma
    repeat = 0;
    neg_z = find(z<=0);
    if isempty(neg_z)==0             % Ensure there are no negative k values
        repeat = 1;
    end
end

k_sp = sparse(u,v,z,Bo,Ao);      % Using SPARSE Matrix
% Adj = spones(k_sp);              % Create Adjacency matrix
k = full(k_sp);                  % Convert from sparse to full matrix

% ************************************************************

% Preallocate memory; Initialize time-dependent sum of molecule numbers
At = zeros(reps,t_steps+1);        Bt = zeros(reps,t_steps+1);        Ct = zeros(reps,t_steps+1);

for n=1:reps

    fprintf(1,'.');
    
    % Arrays for tracking reactant and product agent status
    A = ones(1,Ao);                 B = ones(1,Bo);                 C = zeros(1,Cmax); 
    At(n,1) = sum(A);               Bt(n,1) = sum(B);               Ct(n,1) = sum(C);
    time = zeros(1,t_steps);

    for t=1:t_steps
        
        temp = RandArray(1:size(u,2));      % Randomize order of interacting B-A pairs
        
        for i = 1:size(temp,2)
            
            if B(u(temp(i))) == 1
                   
                    P = 1 - exp(- k(u(temp(i)),v(temp(i))) * dt);     % P_int
%                         P = k(u(temp(i)),v(temp(i))) * dt;            % P_dif

                    if rand < P
                        B(u(temp(i))) = 0;       
                        A(v(temp(i))) = 0;
                        C(u(temp(i))) = 1;   

%                         disp([num2str(u(temp(i))) ',' num2str(v(temp(i)))]);

                    end 
            end
               
        end 

        
    At(n,t+1) = sum(A);         Bt(n,t+1) = sum(B);         Ct(n,t+1) = sum(C);
    time(t+1) = time(t) + dt;
        
    end

end                 % end "for reps" loop

% Population statistics of time trajectories
avgA = mean(At);            avgB = mean(Bt);            avgC = mean(Ct);
sdevA = std(At);            sdevB = std(Bt);            sdevC = std(Ct);

% Solve ODE for 2nd order kinetics
[t_sol,y_sol] = ode45(@(t,y) o2_dif(t,y,ko),0:maxTime/100:maxTime,[Ao ; Bo ; Co]);

%% Plot time trajectories
fig0 = figure('Name','2nd Order Rx Time course','NumberTitle','off'); 
fig0.PaperUnits = 'inches';
fig0.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.

plot(time,avgA,time,avgB);                                              hold on;
% plot(time,avgC);

% *** Plot DE Time Courses ***
plot(t_sol,y_sol(:,1),':c','LineWidth',2);               % scatter(t_sol,y_sol(:,1),'.b');                  
plot(t_sol,y_sol(:,2),':m','LineWidth',2);               % scatter(t_sol,y_sol(:,2),'.g');
% plot(t_sol,y_sol(:,3),':k','LineWidth',2);               % scatter(t_sol,y_sol(:,3),'.r');

p1_devA1 = plot(time,avgA+sdevA,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);
p1_devA0 = plot(time,avgA-sdevA,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);
p1_devB1 = plot(time,avgB+sdevB,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);
p1_devB0 = plot(time,avgB-sdevB,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);
% p1_devC1 = plot(time,avgC+sdevC,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);
% p1_devC0 = plot(time,avgC-sdevC,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);

xlabel('t (sec)');                  ylabel('N(t)');          hold off;
axis([0 maxTime 0 agents+1]);
set(gca,'XMinorTick','on','Box','off');

% leg0 = legend('<N_A(t)>_{sim}','<N_B(t)>_{sim}','DE N_A(t)','DE N_B(t)');
% set(leg0,'FontName','Times New Roman','FontSize',9,...
%     'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');       

% Create textbox
% annotation(fig0,'textbox',[0.02 0.94 0.05 0.07],'String',{'b)'},'FontSize',12,...
%     'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');

% Finish
clear a b z f temp* g i j n x w p* leg* ax* fig*;
toc
