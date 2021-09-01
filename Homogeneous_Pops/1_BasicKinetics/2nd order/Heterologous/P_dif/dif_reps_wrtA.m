%  A + B --> C
% Simulating 2nd order kinetics using my agent-based algorithm.
% The differential rate law for this rx (wrt A) is used
% This is used to calculate the probability of rx for time duration dt.
% The experiment is repeated in order to obtain average trajectories

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

clear; clc;     tic
rng(1);

global k;           % Bimolecular microscopic kinetic constant [units: 1/sec];
k = 0.01;

agents = 50;
maxTime = 10;
dt = 1/100;                               % Fixed time step increment
% dt = 1 / (50*k*agents);                 % Fixed time step increment
% The above equation for dt ensures that the initial probability is <1
% The factor of 50 is arbitrary and assures that P(1)<0.02
reps = 5000;                               % Repeat experiment this # of times
t_steps = ceil(maxTime / dt);

% Preallocate memory; Initialize time-dependent sum of molecule numbers
time = zeros(1,t_steps);
At = zeros(1,t_steps);          Bt = zeros(1,t_steps);          Ct = zeros(1,t_steps);

for n=1:reps

    % Initial number of reactant molecules; assume B is limiting.
    Ao = agents;                    Bo = ceil(0.7*agents);                Cmax = Bo;
    % Arrays for tracking reactant and product agent status
    A = ones(1,Ao);                 B = ones(1,Bo);                 C = zeros(1,Cmax); 
    At(n,1) = sum(A);               Bt(n,1) = sum(B);               Ct(n,1) = sum(C);

    for t=2:t_steps
        temp1 = find(A==1);                 
        P = k * Bt(n,t-1) * dt;         % *** wrt A ***

        for i = 1:size(temp1,2)
            if rand < P
                A(temp1(i)) = 0;                    
                temp2 = find(B==1);                 x = ceil(rand*size(temp2,2));
                B(temp2(x)) = 0;     
                C(temp1(i)) = 1;     
            end
        end 

        time(t) = time(t-1) + dt;
        At(n,t) = sum(A);             Bt(n,t) = sum(B);              Ct(n,t) = sum(C);
        t = t + 1;
    end

end                 % end "for reps" loop

% Population statistics of time trajectories
avgA = mean(At);            avgB = mean(Bt);            avgC = mean(Ct);
sdevA = std(At);            sdevB = std(Bt);            sdevC = std(Ct);

%% Solve ODE for 2nd order kinetics
tmax=time(t-1);
[ty,y_sol] = ode45(@o2_dif,[0:tmax/100:tmax],[At(1) ; Bt(1) ; Ct(1)]);

%% Plot
figure('Name','2nd Order Rx Time course','NumberTitle','off'); 
plot(time,avgA,time,avgB,time,avgC);                                    hold on;

% Plot DE solutions
plot(ty,y_sol(:,1),':b');               % scatter(ty,y_sol(:,1),'.b');                  
plot(ty,y_sol(:,2),':g');               % scatter(ty,y_sol(:,2),'.g');
plot(ty,y_sol(:,3),':r');               % scatter(ty,y_sol(:,3),'.r');

% Plot Std envelopes
p1_devA1 = plot(time,avgA+sdevA,'LineStyle','--','Color',[0.8 0.8 0.8]);
p1_devA0 = plot(time,avgA-sdevA,'LineStyle','--','Color',[0.8 0.8 0.8]);
p1_devB1 = plot(time,avgB+sdevB,'LineStyle','--','Color',[0.8 0.8 0.8]);
p1_devB0 = plot(time,avgB-sdevB,'LineStyle','--','Color',[0.8 0.8 0.8]);
p1_devC1 = plot(time,avgC+sdevC,'LineStyle','--','Color',[0.8 0.8 0.8]);
p1_devC0 = plot(time,avgC-sdevC,'LineStyle','--','Color',[0.8 0.8 0.8]);

xlabel('time');             legend('A','B','C', 'DE A','DE B','DE C');
axis([0 maxTime 0 agents]);

clear p1*;

%% Compare simulation standard deviation with CME...
% obtained from the Chemical Master Equation (see notes in root folder).
% Var = < n^2 > - < n >^2  = < n >           (True in this case; see derivation)
% SDev = sqrt(Var) = sqrt(<n>) 

% Load data from file 'CME_A=??_B=??.mat'    ** First make sure this file exists **
if exist('sigma','var')==0
    cd ../CME_data;
    load(['CME_A=' num2str(Ao) '_B=' num2str(Bo) '.mat'],'sigma');                       
    cd ../P_dif;
    sigma = [0 sigma];
end

% Plot SDev: ABK vs CME
figure2 = figure('Name','SDev Comparison','NumberTitle','off');
set(figure2,'Position',[1 1 500 450]);                          hold on;    
figure2.PaperUnits = 'inches';
figure2.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.

title(['N_{A,i} = ' num2str(Ao) ' , N_{B,i} = ' num2str(Bo)],...
    'FontName','Times New Roman','FontSize',12,'FontWeight','Normal');

plot(time,sdevB,'LineWidth',2);                          
plot(time,sigma,'--','LineWidth',2);
xlabel('t (sec)');           
ylabel('$SDev \Big( <N_B(t)>_{sim} \Big)$','Interpreter','LaTex');                 
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
set(gca,'YTick',[0:0.5:3],'YTickLabel',{'0','0.50','1.00','1.50','2.00','2.50','3.00'});
axis([0 maxTime 0 3]);                                        hold off;

leg2 = legend('ABK','CME');
set(leg2,'Location','SouthEast');
set(leg2,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95]);

% Create textbox
annotation(figure2,'textbox',...
    [0.00771 0.93471 0.05010 0.06429],'String','f)',...
    'FontSize',12,'FontName','Times New Roman',...
    'FitBoxToText','off','EdgeColor','none');

% Find Coefficient of Determination, R^2
SST_s = sum((sdevB - mean(sdevB)).^2);          % Total sum of squares for simulation data (for A)
SSR_s = sum((sdevB - sigma).^2);         % sum of square residuals (sim data vs DE predictions)

Rsq_s = 1 - SSR_s ./ SST_s                      % Definition of R^2

%% Finish
clear temp temp1 temp2 i n x w;
toc

save(['A=' num2str(Ao) '_B=' num2str(Bo) '.mat']);
