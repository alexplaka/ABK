% A --> B + C
% Simulating 1st order (dissociation) kinetics using the Gillespie algorithm.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;
rng(0);

A = ones(1,1000);
lambda = 0.8;
S(1) = 1000;    time(1) = 0;
t = 2;

while S(t-1) > 0
    dt = -log(rand) / ((log(2)/lambda)*S(t-1));
    time(t) = time(t-1) + dt;
    n = ceil(rand*1000);

%     agent = testAgent_r(A,n);   % Choosing only available agents (recursively)
    temp = find(A==1);      x = ceil(rand*size(temp,2));
    agent = temp(x);
    A(agent) = 0;
    S(t) = sum(A);
    
    t = t + 1;
end

% plot(time,S,'r'); hold on;
FitCurve_exp_FixA_g(time',S',1000); hold on;
ezplot(['y=1000*exp(-x*log(2)/' num2str(lambda) ')'],[0 time(t-1) 0 S(1)]);

