%  2A + B --> D
% B is limiting in this scenario. Algorithm is run wrt B.
% Simulating 3rd order kinetics using my agent-based algorithm.
% Assume volume is --- microliters.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

clear; clc;
rng(0);

Avonum = 6.02e+23;      % Avogadro's number
km = 1;      % Bimolecular molar kinetic constant [units: 1 /(M sec)];
V = 10^-21;          % in Liters

global k;       % Bimolecular microscopic kinetic constant [units: 1/sec];
k = km /(Avonum*V)^2;

%% Initial number of reactant molecules; assume B is limiting.
agents = 50;
Ao = agents;      Bo = 0.4*agents;          Dmax = Bo;
% Arrays for tracking reactant and product agent status
A = ones(1,Ao);        B = ones(1,Bo);      D = zeros(1,Dmax); 
% Initialize time-dependent sum of molecule numbers
At(1) = sum(A);        Bt(1) = sum(B);      Dt(1) = 0;

time(1) = 0;
t = 2;
dt = 0.01;                            % Fixed time step increment

while Bt(t-1) > 1                          % Use limiting reactant!
%     dt = 10  / ( At(t-1) * (At(t-1)-1) );       % Variable time step increment

    P = k * At(t-1) * (At(t-1)-1) * dt;
%     P = 1 - exp(-k * At(t-1) * (At(t-1)-1) * dt);
    
    temp = find(B==1);
    for i = 1:size(temp,2)
        if rand < P
            B(temp(i)) = 0;           % Agent B dies
            temp2 = find(A==1);      
            x = ceil(rand(1,2)*size(temp2,2));  % Pick 2 A agents randomly
            A(temp2(x(1))) = 0;        % First agent A dies
            A(temp2(x(2))) = 0;        % Second A agent dies
            D(temp(i)) = 1;            % Agent D is born!
        end
    end 
    time(t) = time(t-1) + dt;
    At(t) = sum(A);         Bt(t) = sum(B);         Dt(t) = sum(D);
    t = t + 1;
end

figure('Name','3rd Order Rx Time course','NumberTitle','off'); 
plot(time,At,'b',time,Bt,'g',time,Dt,'r');             
xlabel('time');                                         hold on;

% Solve ODE for 2nd order kinetics - Plot Time Courses
tmax=time(t-1);
[ty,y_sol] = ode45(@o3_het_dif,[0:tmax/100:tmax],[At(1) ; Bt(1) ; Dt(1)]);
scatter(ty,y_sol(:,1),'xb');                  
scatter(ty,y_sol(:,2),'xg');
scatter(ty,y_sol(:,3),'xr');

legend('A','B','D');        hold off;

clear temp temp2 x i;

% Result:
% - It works for both fixed and variable time step increments. 
