% A + B --> X
% Compute PDF of 2nd order heterologous process from CME solution. 
% Probability expression is Equation 2.99 (p. 68) from reference: 
% Lecca, P. et al. in (book section)
% "Deterministic Versus Stochastic Modelling in Biochemistry and Systems Biology", 
% pp. 35–82 (2013).

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

clear; clc

k = 0.01;

agents = 50;
Ao = agents;                    
Bo = ceil(0.7*agents);      % B is limiting!  
Co = 0;

Z = Ao - Bo;

t = 0.01:0.01:9.99;
tau = k*t;

P = zeros(Bo+1,size(t,2));

for x = 0:Bo
    sum = 0;
    for n = x:Bo
        temp = (-1)^(n-x) * (Z+2*n)*factorial(Z+x+n-1) ...
            / ( factorial(Bo-n)*factorial(n-x)*factorial(Ao+n) ) .* exp(-n*(Z+n).*tau);
        sum = sum + temp;
    end

    P(x+1,:) = factorial(Ao) * factorial(Bo) / (factorial(Z+x)*factorial(x)) * sum;
end

clear temp x n sum;

% Plot PDF of population size of B at time ts
% ts = 15;
% bar(0:Bo,P(:,ts));
% xlabel(['N_B(' num2str(ts) ' sec)']);             ylabel('P_{N_B}');

%% Calculate 1st moment (expectation value) and 2nd moment (variance)
E = zeros(1,size(t,2));                 n_sq = zeros(1,size(t,2));

for a=1:size(t,2)
   for n=0:Bo
      E(1,a) = E(1,a) + n * P(n+1,a); 
      n_sq(1,a) = n_sq(1,a) + n^2 * P(n+1,a); 
   end
end

Var = n_sq - E.^2;                       % calculate variance
sigma = sqrt(Var);                    % calculate standard deviation

% Solve ODE for 2nd order kinetics
[t_sol,y_sol] = ode45(@(t,y) o2_dif(t,y,k),t,[Ao ; Bo ; Co]);

% Plot and compare to DE solution
plot(t,E);                                                                hold on;
plot(t_sol,y_sol(:,2),':g','LineWidth',2);       
plot(t,E + sigma,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);
plot(t,E - sigma,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);
xlabel('t (sec)');                  ylabel('N_B');                        hold off;
axis tight;                         % axis([0 max(t) 0 1]);
set(gca,'XMinorTick','on','Box','off');
legend('CME soln','DE soln')

%% Plot normalized time trajectory
% nE = E/Bo;                          % normalize
% nsigma = sigma / Bo;                % normalize

% plot(nE);                                                                hold on;
% plot(nE + nsigma,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);
% plot(nE - nsigma,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);
% xlabel('t (sec)');                  ylabel('fraction of N_B');           hold off;
% axis tight;                         % axis([0 max(t) 0 1]);
% set(gca,'XMinorTick','on','Box','off');
