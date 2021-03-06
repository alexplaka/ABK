%  A + B --> C
% Simulating 2nd order kinetics using my agent-based algorithm.
% The differential rate law for this rx (wrt B being limiting) is used
% to calculate the probability of rx for time duration dt.
% The experiment is repeated in order to obtain average trajectories

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

clear; clc;     tic

global k;           % Bimolecular microscopic kinetic constant [units: 1/sec];
k = 0.01;

agents = 10;
maxTime = 30;
dt = 1/100;                            % Fixed time step increment
reps = 500;                            % Repeat experiment this # of times
t_steps = ceil(maxTime / dt);

% Preallocate memory; Initialize time-dependent sum of molecule numbers
time = zeros(1,t_steps);
At = zeros(1,t_steps);          Bt = zeros(1,t_steps);          Ct = zeros(1,t_steps);

for n=1:reps

    % Initial number of reactant molecules; assume B is limiting.
    Ao = agents;                    Bo = ceil(0.7*agents);                Cmax = Bo;
    % Arrays for tracking reactant and product agent status
    A = ones(1,Ao);                 B = ones(1,Bo);                 C = zeros(1,Cmax); 
    At(n,1) = sum(A);               Bt(n,1) = sum(B);               Ct(n,1) = sum(C);

    for t=2:t_steps
        temp1 = find(B==1);                 
        P = k * At(n,t-1) * dt;         % *** wrt B ***

        for i = 1:size(temp1,2)
            if rand < P
                B(temp1(i)) = 0;                    
                temp2 = find(A==1);                 x = ceil(rand*size(temp2,2));
                A(temp2(x)) = 0;     
                C(temp1(i)) = 1;     
            end
        end 

        time(t) = time(t-1) + dt;
        At(n,t) = sum(A);             Bt(n,t) = sum(B);              Ct(n,t) = sum(C);
        t = t + 1;
    end

end                 % end "for reps" loop

% Population statistics of time trajectories
avgA = mean(At);            avgB = mean(Bt);            avgC = mean(Ct);
sdevA = std(At);            sdevB = std(Bt);            sdevC = std(Ct);

% Solve ODE for 2nd order kinetics
tmax=time(t-1);
[ty,y_sol] = ode45(@o2_dif,[0:tmax/100:tmax],[At(1) ; Bt(1) ; Ct(1)]);

%% Plot
figure('Name','2nd Order Rx Time course','NumberTitle','off'); 
plot(time,avgA,time,avgB,time,avgC);                                    hold on;
p1_devA1 = plot(time,avgA+sdevA,'LineStyle','--','Color',[0.8 0.8 0.8]);
p1_devA0 = plot(time,avgA-sdevA,'LineStyle','--','Color',[0.8 0.8 0.8]);
p1_devB1 = plot(time,avgB+sdevB,'LineStyle','--','Color',[0.8 0.8 0.8]);
p1_devB0 = plot(time,avgB-sdevB,'LineStyle','--','Color',[0.8 0.8 0.8]);
p1_devC1 = plot(time,avgC+sdevC,'LineStyle','--','Color',[0.8 0.8 0.8]);
p1_devC0 = plot(time,avgC-sdevC,'LineStyle','--','Color',[0.8 0.8 0.8]);

xlabel('time');             legend('Agent A','Agent B','Agent C');
axis([0 maxTime 0 agents]);

%Plot DE solutions
plot(ty,y_sol(:,1),':b');               % scatter(ty,y_sol(:,1),'.b');                  
plot(ty,y_sol(:,2),':g');               % scatter(ty,y_sol(:,2),'.g');
plot(ty,y_sol(:,3),':r');               % scatter(ty,y_sol(:,3),'.r');

clear p1*;
%% Finish
clear temp temp1 temp2 i n x w;
toc

save(['A=' num2str(Ao) '_B=' num2str(Bo) '.mat']);
