% Solution to ODE for a heterogeneous 2nd order chemical reaction:
% A + B --> C
% Note: Concentrations are in numbers of molecules. Microscopic rates and
% constants are considered.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

clear;

Avonum = 6.02e+23;          % Avogadro's number        
V = 10^-21;                 % Volume in Liters
km = 0.01;       % Bimolecular molar kinetic constant [units: 1 /(M sec)];
global k;       
k = km /(Avonum*V); % Bimolecular microscopic kinetic constant [units: 1/sec];

A(1) = 1000;        % Initial particle number
B(1) = 700;         % Initial particle number  
C(1) = 0;           % Initial particle number

% Am = A(1)/Avonum/V;     Bm = B(1)/Avonum/V;
tmax = 400;

[t y_sol] = ode45(@o2_dif,[0:tmax/100:tmax],[A(1) ; B(1) ; C(1)]);
figure('Name','2nd Order Rx Time course','NumberTitle','off'); 
scatter(t,y_sol(:,1),'.b');                  hold on;
scatter(t,y_sol(:,2),'.g');
scatter(t,y_sol(:,3),'.r');
xlabel('time');             legend('A','B','C');

% Result: 
% - For such small initial numbers of molecules, the volume has 
% to be small enough for the reaction to progress fast enough within the 
% enforced time frame (tmax = 400).
% - It works (given the adjustment to a microscopic rate constant)!