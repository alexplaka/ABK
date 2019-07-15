% Probability of discrete rx events for 1st order process.
% Comparison between exact P and its approximation as a Poisson distribution.
% See Report: sec 2.2.5 "From Individuals to Populations"

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

clear;                  % clc;

Na = 10;                % number of particles
k= 0.1;                 % 1st order rate constant
dt = 0.01;              % time interval

r = k * Na * dt;
mu_max = 4;

for mu=1:mu_max                  % mu: reaction events occurring in dt
   Poisson(mu) = r^mu / factorial(mu) * exp(-r);        % Poisson probability
   P(mu) = nchoosek(Na,mu) * (1 - exp(-k*dt))^mu * (exp(-k*dt))^(Na-mu);
   % Full probability expression
end

% *** Note: mu = 0 gives identical values for P and Poisson ***

% plot(1:mu_max,Poisson,'b');                         hold on;
% plot(1:mu_max,P,'.r');
% xlabel('\mu');                      ylabel('P');

percDiff = abs(Poisson-P)./P * 100;

disp(['Na = ' num2str(Na) ' , k = ' num2str(k) ' , dt = ' ...
    num2str(dt) ' , mu = 1:' num2str(mu_max)]);

disp(['% diff = ' num2str(percDiff)]);

% Result:
% The Poisson distribution is a reasonable approximation only when Na>100 and small dt.
% For smaller populations, dt has to be even smaller such that only 0 or 1 reaction events
% are possible within dt.