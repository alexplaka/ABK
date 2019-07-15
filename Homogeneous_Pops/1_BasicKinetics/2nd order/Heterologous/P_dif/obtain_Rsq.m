clear;

Ao = 10;                  Bo = 7;                  Co = 0;

cd('wrtB_data');
filename = ['A=' num2str(Ao) '_B=' num2str(Bo) '.mat'];
load(filename);             cd ..;

[t_sol, y_sol] = ode45(@o2_dif,time,[Ao;Bo;Co]);

%% Find Coefficient of Determination, R^2
% This works for any curve fit, linear or nonlinear.
SST = sum((avgB' - mean(avgB)).^2);       % Total sum of squares for simulation data (for B)

SSR = sum((avgB' - y_sol(:,2)).^2);       % sum of square residuals (sim data vs DE predictions)
Rsq = 1 - SSR/SST                         % Definition of R^2
R = sqrt(Rsq);                            % Correlation Coefficient R

%% Finish
cd('wrtB_data');
save(filename,'Rsq','-append');
cd ..;
