% Compute PDF of 2nd order rx from CME solution (HOMOGENEOUS POPULATION). 
% Probability expression is Equation 2.99 (p. 68) from reference: (book section)
% Lecca, P. et al. in Deterministic Versus Stochastic Modelling ...
% in Biochemistry and Systems Biology, pp. 35–82 (2013).

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

k = 0.005;
agents = 10;
Ao = agents;                    Bo = ceil(0.7*agents);      % B is limiting!  
Co = 0;
Z = Ao - Bo;
t = 1:200;
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

clear x n sum temp;
%% Plot PDF of population size of B at time ts
ts = 200;
bar(0:Bo,P(:,ts));
xlabel(['N_B(' num2str(ts) ' sec)']);             ylabel('P_{N_B}');

%% Calculate 1st moment (expectation value) and 2nd moment (variance)
E = zeros(1,size(t,2));                 V = zeros(1,size(t,2));

for a=1:size(t,2)
   for n=0:Bo
      E(1,a) = E(1,a) + n * P(n+1,a); 
      V(1,a) = V(1,a) + n^2 * P(n+1,a); 
   end
end

Var = V - E.^2;
sigma = sqrt(Var);                    % calculate standard deviation

clear a n;
%% Plot 
plot(E);                                                                hold on;
plot(E + sigma,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);
plot(E - sigma,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);
xlabel('t (sec)');                  ylabel('N_B');                      hold off;
% axis tight;                         
axis([0 max(t) 0 10]);
set(gca,'XMinorTick','on','Box','off');
legend('CME soln');

%% Plot normalized time trajectory
nE = E/Bo;                          % normalize
nsigma = sigma / Bo;                % normalize

% plot(nE);                                                                hold on;
% plot(nE + nsigma,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);
% plot(nE - nsigma,'LineStyle','--','Color',[0.8 0.8 0.8],'LineWidth',2);
% xlabel('t (sec)');                  ylabel('fraction of N_B');           hold off;
% axis tight;                         % axis([0 max(t) 0 1]);
% set(gca,'XMinorTick','on','Box','off');
