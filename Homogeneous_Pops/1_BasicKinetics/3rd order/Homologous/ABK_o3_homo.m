%  3A --> D
% Simulating 3rd order kinetics using my agent-based algorithm.
% The integrated rate law for this rx is:
% A/Ao = sqrt( 1 / (2 * 3 * k * Ao^2 * t + 1))
% This is used to calculate the probability of rx for time duration dt.
% Notice that the stoichiometric coeeficient, 3, is omitted from the
% the expression for P, since P is evaluated for each molecule
% separately from the others, and the correct number of molecules is
% marked as reacted when rand < P.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear; clc;                             tic
rng(0);

global k dif_canon;       % Bimolecular microscopic kinetic constant [units: 1/sec];

dif_canon = 0;

% Avonum = 6.02e+23;      % Avogadro's number
% km = 1;      % Bimolecular molar kinetic constant [units: 1 /(M sec)];
% V = 10^-21;          % Volume in Liters
% k = km /(Avonum*V)^2;

k = 0.001;
%% Initial number of reactant molecules
agents = 100;
Ao = agents;                                Dmax = Ao;
% Arrays for tracking reactant and product agent status
A = ones(1,Ao);                             D = zeros(1,Dmax); 
% Initialize time-dependent sum of molecule numbers
At(1) = sum(A);                             Dt(1) = 0;

time(1) = 0;
t = 2;
dt = 1/100;                            % Fixed time step increment
% dt = 1 / (20*k*Ao^2);                   % Fixed time step increment
% The above equation for dt ensures that the initial probability is <1
% The factor of 20 is arbitrary and assures that P(1)<0.05

while At(t-1) > agents/10                     
%     dt = 10  / ( (At(t-1)-1) * (At(t-1)-2) );       % Variable time step increment

%     P(t-1) = k * (At(t-1)-1) * (At(t-1)-2) * dt;     % Differential form
    % Integrated form of P:
    P(t-1) = 1 - sqrt(1 / (2 * k * (At(t-1)-1) * (At(t-1)-2) * dt + 1));
    
    temp = find(A==1);
    for i = 1:size(temp,2)
        if A(temp(i))==1 && rand < P(t-1)
        % First IF condition ensures that an A agent has not already "died" in this time step
            A(temp(i)) = 0;           % First agent A dies
            temp2 = find(A==1);      
            x = ceil(rand(1,2)*size(temp2,2));  % Pick 2 A agents randomly
            A(temp2(x(1))) = 0;        % Second agent A dies
            A(temp2(x(2))) = 0;        % Third agent A dies
            D(temp(i)) = 1;            % Agent D is born!
        end
    end 
    time(t) = time(t-1) + dt;
    At(t) = sum(A);                     Dt(t) = sum(D);
    t = t + 1;
end

figure('Name','3rd Order Rx Time course','NumberTitle','off'); 
plot(time,At,'b',time,Dt,'r');             
xlabel('time');                                         hold on;

% Solve ODE for 3rd order kinetics - Plot Time Courses
tmax=time(t-1);
[ty,y_sol] = ode45(@o3_homo_dif,[0:tmax/100:tmax],[At(1) ; Dt(1)]);
scatter(ty,y_sol(:,1),'xb');                  
scatter(ty,y_sol(:,2),'xr');

legend('A','D');                            hold off;

clear temp temp2 x i;
toc

% Result: 
% - It works for both fixed and variable time step increments. 
% - Works for both the differential and integrated rate law forms
% of the transition probability.
% - Simulation is a little faster when using the integrated form of P
% (surprisingly, compared to previous observations for 1st order rxs)
% by approx. 15%. 
