%  A + B --> C
% Simulating 2nd order kinetics using ABK.
% Heterogeneous Population!
% Here, I evaluate each possible interaction wrt B.
% (ie, I draw a random number for each possible B-A pair).

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear; clc;
rng(0);

agents = 100;
dt = 0.01;                            % Fixed time step increment

% Initial number of reactant molecules; assume B is limiting.
Ao = agents;            Bo = 0.7*agents;            Cmax = Bo;

% Set up k-matrix wrt B (size: Bo x Ao)
k1 = 0.001*ones(Bo,Ao/2);                   k2 = 0.002*ones(Bo,Ao/2);
k = [k1, k2];
ko = mean(mean(k));

% Arrays for tracking reactant and product agent status
A = ones(1,Ao);        B = ones(1,Bo);          C = zeros(1,Cmax); 
% Preallocate memory and Initialize time-dependent sum of molecule numbers
At(1) = sum(A);        Bt(1) = sum(B);          Ct(1) = 0;

time(1) = 0;
t = 1;

while Bt(t) > 1                         % Use limiting reactant!

    for i = 1:size(B,2)
        if B(i) == 1                    % if this B agent is alive
            for j = 1:size(A,2)         % then see if it reacts with an A agent
                if A(j) == 1            % if this A agent is alive
                    P = 1 - exp(- k(i,j) * dt);
                    if rand < P
                        B(i) = 0;       
                        A(j) = 0;
                        C(i) = 1;   
                    end
                end
            end
        end
    end 
    
    At(t+1) = sum(A);         Bt(t+1) = sum(B);         Ct(t+1) = sum(C);
    
    time(t+1) = time(t) + dt;
    t = t + 1;
end

%%
figure('Name','2nd Order Rx Time course','NumberTitle','off'); 
plot(time,At,time,Bt,time,Ct);             hold on;
xlabel('time');             legend('Agent A','Agent B','Agent C');

% Solve ODE for 2nd order kinetics - Plot Time Courses
tmax=time(t-1);
[ty,y_sol] = ode45(@(t,y) o2_dif(t,y,ko),[0:tmax/100:tmax],[At(1) ; Bt(1) ; Ct(1)]);
scatter(ty,y_sol(:,1),'.b');                  
scatter(ty,y_sol(:,2),'.g');
scatter(ty,y_sol(:,3),'.r');

clear temp x i;

% Result: It works! 