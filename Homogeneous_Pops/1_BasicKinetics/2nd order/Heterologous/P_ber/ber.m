%  A + B --> C
% Simulating 2nd order kinetics using my agent-based algorithm.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

clear; clc;

Avonum = 6.02e+23;      % Avogadro's number
km = 0.01;      % Bimolecular molar kinetic constant [units: 1 /(M sec)];
V = 10^-21;          % Volume in Liters
agents = 10;
dt = 0.01;                            % Fixed time step increment

global k;       % Bimolecular microscopic kinetic constant [units: 1/sec];
k = km /(Avonum*V);

% Initial number of reactant molecules; assume B is limiting.
Ao = agents;      Bo = 0.7*agents;       Cmax = Bo;
% Arrays for tracking reactant and product agent status
A = ones(1,Ao);        B = ones(1,Bo);          C = zeros(1,Cmax); 
% Initialize time-dependent sum of molecule numbers
At(1) = sum(A);        Bt(1) = sum(B);          Ct(1) = 0;

time(1) = 0;
t = 2;

while Bt(t-1) > 1                     % Use limiting reactant!
%     dt = 10  / At(t-1);             % Variable time step increment
%     dt = 1/( k*At(t-1)*Bt(t-1) );   % Variable time step increment (ver. 2)
%     dt = exprnd(1/( k*At(t-1)*Bt(t-1) ) )   % Exponentially-distributed time step increment

    P(t-1) = 1 - exp(-k * At(t-1) * dt);

    for i = 1:size(B,2)
        if B(i) == 1 && rand < P(t-1)
            B(i) = 0;       
            temp = find(A==1);      x = ceil(rand*size(temp,2));
            A(temp(x)) = 0;
            C(i) = 1;   
        end
    end 
    time(t) = time(t-1) + dt;
    At(t) = sum(A);     Bt(t) = sum(B);
    Ct(t) = sum(C);
    t = t + 1;
end

figure('Name','2nd Order Rx Time course','NumberTitle','off'); 
plot(time,At,time,Bt,time,Ct);             hold on;
xlabel('time');             legend('Agent A','Agent B','Agent C');

% Solve ODE for 2nd order kinetics - Plot Time Courses
tmax=time(t-1);
[ty,y_sol] = ode45(@o2_dif,[0:tmax/100:tmax],[At(1) ; Bt(1) ; Ct(1)]);
scatter(ty,y_sol(:,1),'.b');                  
scatter(ty,y_sol(:,2),'.g');
scatter(ty,y_sol(:,3),'.r');

clear temp x i;

% Result: It works! 
% - It works for both fixed and all variable time step increment algorithms. 
% - Simulation is much slower, overall, for fixed time steps (as expected).
