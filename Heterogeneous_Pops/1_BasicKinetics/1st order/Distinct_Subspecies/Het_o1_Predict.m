function [N_avg, N_sdev] = Het_o1_Predict(time_array, Ao, k_array)        
% This function calculates the first and second (central) moments of
% the probability destribution describing the population size of a 
% species undergoing a 1st order transition, A --> X. It takes into 
% account the time-independent PHM of a heterogeneous population, 
% thus allowing the prediction of the average time trajectory and 
% its standard deviation, both of which cannot be computed using
% standard deterministic or stochastic methods.
% Inputs: 
% - time_array: horizontal array of sampled time values.
% - Ao: Initial population size of species A.
% - k_array: horizontal (time-independent) PHM of species A.

% Note 1: I have used the terminology and naming convention of
% "alive" and "dead" A agents, who are still in existence or 
% have already decayed/transitioned, respectively.
% Note 2: This approach is computationally feasible for small 
% population sizes of <20 agents because of the cost of 
% enumerating all possible combinations of alive and dead
% agents (using the Matlab-provided "combnk" function; 
% Documentation suggests not using this function for arrays with
% more than 15 elements).

% I have taken a very careful approach to the calculations that
% are the purpose of this function. Some of the code can be
% condensed, but I have kept each step of the process separate
% for readibility, clarity, and debugging purposes.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

time = time_array;                  
agents = Ao;                        % Initial population size
k = k_array;                        % Time-independent PHM
%% Possible combinations of alive and dead agents
C_alive = cell(1,agents);                               % C_alive_none = [];

for x=1:agents
    C_alive{x}= combnk(1:agents,x);
%     disp(['combnk: alive=' num2str(x)]);          % For debugging purposes
end
clear x;

% Possible combinations of dead agents ("complement" of C)
C_dead = cell(1,agents);                                % C_dead_all = 1:agents;

for a=1:size(C_alive,2)
%     disp(['combnk: alive=' num2str(a)]);
    for b=1:size(C_alive{1,a},1)
        all = 1:agents;
        all(C_alive{1,a}(b,:)) = [];
        C_dead{1,a}(b,:) = all;
    end
%     disp(['combnk: alive=' num2str(a) ' --> dead ; DONE!']);
end

% disp('Comb, dead: done');
clear a b c;
% Note: indexing of C_dead is wrt # of alive agents. 
% e.g., C{1} lists the number of possible combinations of 9 dead agents.
%% Calculate probabilities for all of these ALIVE (and DEAD) agent combinations
% P [alive] = exp(-k_i * t)
P_alive = cell(1,agents,size(time,2));                  % P_alive_none = [];
% Using a singleton dimension for consistency's sake.

for t=1:size(time,2)
    for d = 1:agents
        for f=1:size(C_alive{1,d},1)
            for g=1:size(C_alive{1,d},2)
                P_alive{1,d,t}(f,g) = exp(-k(C_alive{1,d}(f,g)) * time(t));
            end
        end
    end
end
clear t d f g;

% Now calculate probabilities for all of these DEAD agent combinations
% P [dead] = 1 - exp(-k_i * t)
P_dead = cell(1,agents,size(time,2));                  P_dead_all = zeros(1,agents,size(time,2));
% Using a singleton dimension for consistency's sake.

for t=1:size(time,2)
    for h = 1:agents
        
        P_dead_all(1,h,t) = 1 - exp(-k(h) * time(t));
        
        for i=1:size(C_dead{1,h},1)
            for j=1:size(C_dead{1,h},2)
                P_dead{1,h,t}(i,j) = 1 - exp(-k(C_dead{1,h}(i,j)) * time(t));
            end
        end
    end
end
% disp('Calc P: done');
clear t h i j;
%% Calculate products of ALIVE (and DEAD) probabilities within each row of each cell entry. 
P_a = cell(1,agents,size(time,2));                      % P_a_none = zeros(1,size(time,2));

for t=1:size(time,2)
    for m = 1:agents
        for n=1:size(C_alive{1,m},1)
            P_a{1,m,t}(n,1) = prod(P_alive{1,m,t}(n,:));
        end
    end
end
clear t m n;

% Now calculate products of DEAD probabilities within each row of each cell entry. 
P_d = cell(1,agents,size(time,2));                      P_d_all = zeros(1,size(time,2));

for t=1:size(time,2)
    for p = 1:agents-1                   % b/c C_dead{1,end} = [];
        for q=1:size(C_dead{1,p},1)
            P_d{1,p,t}(q,1) = prod(P_dead{1,p,t}(q,:));
        end
    end
    P_d{1,agents,t} = [];
    P_d_all(1,t) = prod(P_dead_all(1,:,t));
end
% disp('Products of P: done');
clear t p q;

%% Multiply probabilities for each combination of alive and dead agents 
% --> Total probability for each combination.
P_ad = cell(1,agents,size(time,2));                    

for t=1:size(time,2)
    for r = 1:agents-1                  % b/c P_d{1,agents,t} = [];
        for s=1:size(P_a{1,r,t},1)
            P_ad{1,r,t}(s,1) = P_a{1,r,t}(s,1) * P_d{1,r,t}(s,1);
        end
    end
    P_ad{1,agents,t} = P_a{1,agents,t};
end
% disp('Calc P_tot for each comb: done');
clear t r s;

%% Sum probabilities all possible combinations for each value of alive agents (1:agents)
P = zeros(agents,size(time,2));                        
% P_d_all stores the probabilities of zero (0) alive agents as a function of time.

for t=1:size(time,2)
    for u = 1:agents
        P(u,t) = sum(P_ad{1,u,t});
    end
end
% disp('Sum P_tot of all combs: done')
clear t u;

%% Obtain expected time trajactory
N_avg = zeros(1,size(time,2));          % Expected (avg) number of alive agents (function of time).
N_m2 = zeros(1,size(time,2));           % Raw second moment of alive agents (function of time).

for t=1:size(time,2)
    N_avg(t) = 0 * P_d_all(t);          % written in full form for clarity
    N_m2(t) = 0^2 * P_d_all(t); 
    
    for v=1:agents
       N_avg(t) = N_avg(t) + v * P(v,t); 
       N_m2(t) = N_m2(t) + v^2 * P(v,t); 
    end
end

N_var = N_m2 - N_avg.^2;
N_sdev = sqrt(N_var);

% disp('Calc moments: done');
