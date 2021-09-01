%  A --> B
% Simulating 1st order kinetics. Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

clear;       tic;
rng(0);
    
totalTime = 10;             % Simulation time (sec)
% dt = 1/20;                % Constant (fixed) time step increment (sec)

if exist('dt') == 1         
    t_steps = totalTime / dt;
else
    t_steps = 1200;             % empirically chosen value (for lambda = 0.8)
end

agents = 1000;
A = ones(t_steps,agents);      
S = zeros(1,t_steps);           S(1) = agents;
lambda = 0.8;                   k = log(2)/lambda;
time = zeros(1,t_steps);
P = zeros(1,t_steps);

for t = 2:t_steps
    dt = 1 / (k*S(t-1));            % Variable time step increment
%     dt = exprnd(1/(k*S(t-1)));   % Exponentially-distributed time step increment
    time(t) = time(t-1) + dt;      
    for i = 1:size(A,2)
        if A(t,i) == 1                      % if agent is still alive
            P(t) = 1 - 2^(-dt/lambda);      % P from integrated rate law
%             P(t) = k*dt;                  % P from differential rate law
            if rand < P(t)                  % **check probability condition**
                A(t,i) = 0;                 % agent dies
                for w = t+1:t_steps
                    A(w,i) = 0;             % mark agent dead for rest of time
                end
            end
        end
    end
    S(t) = sum(A(t,:));
end

k_est = FitCurve_exp_FixA_g(time',S',agents);
hold on;
ezplot(['y=' num2str(agents) '*2^(-x/' num2str(lambda) ')'], ...
    [0 time(t_steps) 0 agents]);
xlabel('t');            ylabel('A');

clear t w i;
toc
