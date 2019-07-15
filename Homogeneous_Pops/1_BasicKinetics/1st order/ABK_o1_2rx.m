% Simulating 1st order kinetics. Considering individual agents.
% Here, I consider the case where a single reactant species A undergoes two
% reactions, each with the indicated rate constant.
% A --> B       k1 = 0.80 /sec 
% A --> C       k2 = 0.29 /sec
% (This is meant to simulate the reactions a CaM-bound CaMKII subunit can
% undergo: CaM dissociation and autophosporylation (considered to be
% pseudo-first order for now, and used merely as an example).

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

clear;                    tic;

% Define Parameters
t_steps = 1200;              agents = 1000;
A = ones(t_steps,agents);    B = zeros(t_steps,agents);   C = zeros(t_steps,agents);
Asum = zeros(1,t_steps);     Bsum = zeros(1,t_steps);     Csum = zeros(1,t_steps);       
Asum(1) = agents;            Bsum(1) = 0;                 Csum(1) = 0;
k1 = 0.80;                   k2 = 0.29;      
lambda1 = log(2)/k1;         lambda2 = log(2)/k2;  

% dt = 1/100;               % Constant time step increment
time(1) = 0;

for t = 2:t_steps
    dt = 1 / ((k1+k2)*Asum(t-1));            % Variable time step increment
%     dt = exprnd(1/( (k1+k2)*Asum(t-1) ) );   % Exponentially-distributed time step increment
    time(t) = time(t-1) + dt;
    for i = 1:size(A,2)
        if A(t,i) == 1                 % if agent is still alive
%             P(t) = 1 - 2^(-dt*(1/lambda1 + 1/lambda2));
%             P1(t) = P(t) * lambda2/(lambda1+lambda2);
%             P2(t) = P(t) * lambda1/(lambda1+lambda2);
            P(t) = 1 - exp(-dt*(k1+k2));
            P1(t) = P(t) * k1/(k1+k2);
            P2(t) = P(t) * k2/(k1+k2);
            r = rand;
            if r < P1(t)   % ****check probability condition for Rx 1
                A(t,i) = 0;            % agent death
                B(t,i) = 1;            % agent birth
                for w = t+1:t_steps    A(w,i) = 0;      end % mark agent dead for rest of time
                for w = t+1:t_steps    B(w,i) = 1;      end % mark agent alive
            elseif r >= P1(t) && r < P1(t)+P2(t)  % ****check probability condition for Rx 2
                A(t,i) = 0;            % agent death
                C(t,i) = 1;            % agent birth
                for w = t+1:t_steps    A(w,i) = 0;      end % mark agent dead for rest of time
                for w = t+1:t_steps    C(w,i) = 1;      end % mark agent alive
            end
        end
    end
    Asum(t) = sum(A(t,:));
    Bsum(t) = sum(B(t,:));
    Csum(t) = sum(C(t,:));
end

lambda_theory = log(2)/(k1+k2)
%%
% plot(S);
FitCurve_exp_FixA_g(time',Asum',agents);
hold on;
ezplot(['y=' num2str(agents) '*exp(-' num2str(k1+k2) '*x)'],[0 5 0 agents]);

toc

% Result: 
% 1st order kinetics for the disappearance of A is observed with 
% k = k1 + k2, using both constant and variable time steps.
% Note that overall lambda = log(2)/(k1+k2) 
% In this case lambda = 0.6359 (theoretically).