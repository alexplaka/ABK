%  A + B --> C
% Simulating 2nd order kinetics using the Gillespie algorithm
% for a heterologous 2nd order process.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear;
rng(0);

Avonum = 6.02e+23;      % Avogadro's number
km = 0.01;      % Bimolecular molar kinetic constant [units: 1 /(M sec)];
V = 10^-21;          % in Liters

global k;       % Bimolecular microscopic kinetic constant [units: 1/sec];
k = km /(Avonum*V);

A(1) = 1000;                                        % Initial particle number  
B(1) = 700;                                         % Initial particle number
C(1) = 0;                                           % Initial particle number
%  Am = A(1)/Avonum/V;     Bm = B(1)/Avonum/V;      % Molar concentrations

time(1) = 0;
t = 2;

while B(t-1) > 0                                    % Use limiting reactant!
    a(t) = k*A(t-1)*B(t-1);                         % Propensity function
    dt = -log(rand) / a(t);
    time(t) = time(t-1) + dt;
    A(t) = A(t-1) - 1;
    B(t) = B(t-1) - 1;
    C(t) = C(t-1) + 1;    
    t = t + 1;
end

figure('Name','2nd Order Rx Time course','NumberTitle','off'); 
plot(time,A,time,B,time,C);                 hold on;
xlabel('time');             legend('Gil A','Gil B','Gil C');

% Solve ODE for 2nd order kinetics - Plot Time Courses
tmax=time(t-1);
[t y_sol] = ode45(@o2_het_dif,[0:tmax/100:tmax],[A(1) ; B(1) ; C(1)]);
scatter(t,y_sol(:,1),'.b');                  
scatter(t,y_sol(:,2),'.g');
scatter(t,y_sol(:,3),'.r');


