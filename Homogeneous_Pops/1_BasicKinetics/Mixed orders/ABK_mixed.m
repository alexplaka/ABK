%  A     --> C              Rate constant k1  (1st order)
%  A + B --> D              Rate constant k2  (2nd order)
% Simulating concurrent 1st and 2nd order kinetics using my agent-based
% algorithm. Assume volume is --- microliters.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear; clc;  tic
rng(0);

Avonum = 6.02e+23;      % Avogadro's number
V = 10^-21;          % in Liters

k2m = 1;      % Bimolecular **molar** kinetic constant [units: 1 /(M sec)];
% k1: 1st order rate constant [units: 1/sec]
% k2: bimolecular **microscopic** kinetic constant [units: 1/sec]
global k1 k2;
k1 = 0.8;
k2 = k2m /(Avonum*V);

% Initial number of reactant molecules; Assume B is limiting for 2nd order rx.
Ao = 100;      Bo = 70;       
% Arrays for tracking reactant agent status
A = ones(1,Ao);                 B = ones(1,Bo);          
% Initialize time-dependent sum of molecule numbers
At(1) = sum(A);        Bt(1) = sum(B);          
Ct(1) = 0;             Dt(1) = 0;

time(1) = 0;
t = 2;
dt = 0.01;                                      % Fixed time step increment

while Bt(t-1) > 0 && At(t-1) > 0
%     dt = 10 / At(t-1);                        % Variable time step incr. 1
%     dt = 1/(k1*At(t-1) + k2*At(t-1)*Bt(t-1)); % Variable time step incr. 2
%     dt = exprnd(1/(k1*At(t-1) + k2*At(t-1)*Bt(t-1)));   % Exponentially-distributed time step increment
    P1 = k1 * dt;             % Prob of 1st order rx A --> C
%     P1(t-1) = 1 - exp(-k1 * dt);    % Prob of 1st order rx A --> C
    P2(t-1) = k2 * Bt(t-1) * dt;   % Prob of 2nd order rx A + B --> D
%     Ptot(t-1) = P1(t-1) + P2(t-1);
    deltaC = 0;             deltaD = 0;  
    for i = 1:size(A,2)
        if A(i) == 1      
            r = rand;
            if r < P1          % 1st order rx
                A(i) = 0;           deltaC = deltaC + 1;
            elseif r >= P1 && r < P1+P2(t-1)  % 2nd order rx
                A(i) = 0;           
                temp = find(B==1);      x = ceil(rand*size(temp,2));
                B(temp(x)) = 0;
                deltaD = deltaD + 1;
            end 
        end
    end 
    At(t) = sum(A);                 Bt(t) = sum(B);
    Ct(t) = Ct(t-1) + deltaC;       Dt(t) = Dt(t-1) + deltaD;
    time(t) = time(t-1) + dt;
    t = t + 1;
end

%% Plot trajectories
figure('Name','Mixed Rx Time Course','NumberTitle','off'); 
plot(time,At,time,Bt,time,Ct,time,Dt);            hold on;
xlabel('time');       

% Solve ODE for mixed order kinetics - Plot Time Courses
tmax=time(t-1);
[ty,y_sol] = ode45(@o12_mixed_dif,[0:tmax/100:tmax],[At(1);Bt(1);Ct(1);Dt(1)]);
scatter(ty,y_sol(:,1),'.b');         
scatter(ty,y_sol(:,2),'.g');
scatter(ty,y_sol(:,3),'.r');
scatter(ty,y_sol(:,4),'.c');

legend('A','B','C','D','DE A','DE B','DE C','DE D','Location','Best');

%%
clear temp x deltaC deltaD r i;
toc

% Result: 
% - Works for fixed, variable, and exponentially-distributed time step
% increments. (latter has the smallest computing time).
% - Note that calculated probabilities work best for small dt. For variable
% time step increments some dt's will be large (towards the end of the
% simulation). It may be best to use fixed time step increments while
% calculating such probabilities, despite its computational cost.