% Solution to ODE for a heterogeneous 2nd order chemical reaction:
% A + B --> C
% Note: Concentrations are molarities.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

clear;

Avonum = 6.02e+23;          % Avogadro's number        
V = 10^-21;                 % Volume in Liters

global k;       % Bimolecular molar kinetic constant [units: 1 /(M sec)];
k = 0.01;

Ap = 1000;          % Initial particle number 
Bp = 700;           % Initial particle number
Cp = 0;             % Initial particle number
A(1) = Ap/(Avonum*V);     B(1) = Bp/(Avonum*V);     C(1) = Cp/(Avonum*V);
tmax = 400;

[t y_sol] = ode45(@o2_dif,[0:tmax/100:tmax],[A(1) ; B(1) ; C(1)]);
figure('Name','2nd Order Rx Time course','NumberTitle','off'); 
scatter(t,y_sol(:,1),'.b');                  hold on;
scatter(t,y_sol(:,2),'.g');
scatter(t,y_sol(:,3),'.r');
xlabel('time');             legend('A','B','C');