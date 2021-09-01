%  A + A --> C
% Simulating 2nd order kinetics using the Gillespie algorithm
% for a homologous 2nd order process.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;
rng(0);

% Avonum = 6.02e+23;      % Avogadro's number
% km = 0.01;      % Bimolecular molar kinetic constant [units: 1 /(M sec)];
% V = 10^-21;          % in Liters

global k;       % Bimolecular microscopic kinetic constant [units: 1/sec];

% k = km /(Avonum*V);
k = 0.01;

A(1) = 100;                                        % Initial particle number  
C(1) = 0;                                           % Initial particle number

time(1) = 0;
t = 2;

while A(t-1) > 0                                    
    a(t) = k/2 * A(t-1) * (A(t-1)-1);                  % Propensity function
    dt = -log(rand) / a(t);
    time(t) = time(t-1) + dt;
    A(t) = A(t-1) - 2;
    C(t) = C(t-1) + 1;    
    t = t + 1;
end

figure('Name','2nd Order Rx Time course','NumberTitle','off'); 
plot(time,A);                                        hold on;
% plot(time,C);
xlabel('time');             
legend('Gil A');

% Solve ODE for 2nd order kinetics - Plot Time Courses
tmax=time(t-1);
[t y_sol] = ode45(@o2_hom_dif,[0:tmax/100:tmax],[A(1) ; C(1)]);
scatter(t,y_sol(:,1),'.b');                  
% scatter(t,y_sol(:,2),'.r');

