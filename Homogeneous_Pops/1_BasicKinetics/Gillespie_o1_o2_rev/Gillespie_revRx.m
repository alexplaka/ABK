% This is an implementation of the Gillespie algorithm for the stochastic
% simulation of a simple reversible reaction A <--> B.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

clear;
rng(0);

k1 = 0.1;                   % Define rate constant for forward reaction
k2 = 0.2;                   % Define rate constant for reverse reaction
A(1) = 100;   B(1) = 0;     % Initialize numbers of molecules
t = 0;    time(1) = 0;      % Initialize time
t_max = 100;                % Maximum time for simulation to run
n = 1;

while t < t_max
    v1 = k1 * A(n);
    v2 = k2 * B(n);
    h = v1 + v2;
    nv1 = v1 / h;
    nv2 = v2 / h;

    p = rand(1,2);          % Two uniformly distributed random numbers
    
    tau = - log(p(1)) / h;  
    t = t + tau;        	% When next reaction will occur
    
    if p(2) < nv1           % Which reaction will occur
        A(n+1) = A(n) - 1;      B(n+1) = B(n) + 1;
    else
        A(n+1) = A(n) + 1;      B(n+1) = B(n) - 1;
    end
    
    n = n + 1;
    time(n) = t;
end

plot(time,A,time,B);
title('Stochastic simulation using Gillespie algorithm: A <--> B');
xlabel('time'); ylabel('# of molecules');
legend('A','B');

clear v1 v2 h nv1 nv2 p;