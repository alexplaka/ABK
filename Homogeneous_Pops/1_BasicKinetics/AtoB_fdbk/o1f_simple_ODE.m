% Solution to ODE for a 1st order chemical reaction with feedback:
% A --> B  ,  where B influences the rate of the reaction.

% Note: Concentrations are in numbers of molecules. Microscopic rates and
% constants are considered.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear;

global ko alpha K n;

lambda = 0.8;               % half-life of first order process     
ko = log(2)/lambda;         % basal rate: unimolecular kinetic constant [units: 1 / sec]
% Feedback parameters
alpha = 0;                  % degree of activation: >1 activator, <1 repressor, =1 no regulation
n = 1;                      % Hill Coefficient: measure of cooperativity
K =  50;                    % Half-maximal activation/repression

A(1) = 100;         % Initial particle number
B(1) = 0;           % Initial particle number  

tmax = 100;

[t y_sol] = ode45(@o1f_simple_dif,[0:tmax/100:tmax],[A(1) ; B(1)]);
figure('Name','1st Order Rx w/ Feedback Time course','NumberTitle','off'); 
scatter(t,y_sol(:,1),'.b');                  hold on;
scatter(t,y_sol(:,2),'.g');

xlabel('time');             legend('A','B');
