%  A + B --> C
% Simulating 2nd order kinetics using my agent-based algorithm.
% The differential rate law for this rx (wrt B being limiting) is used
% This is used to calculate the probability of rx for time duration dt.
% The experiment is repeated in order to obtain average trajectories

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

clear; clc;     tic

global k;           % Bimolecular microscopic kinetic constant [units: 1/sec];
k = 0.01;

agents = 10;
maxTime = 30;
dt = 1/100;                            % Fixed time step increment
t_steps = ceil(maxTime / dt);

reps = 2000;                              % Repeat experiment this # of times

% Initial number of reactant molecules; assume B is limiting.
Ao = agents;                    Bo = ceil(0.7*agents);                Co = 0;
Cmax = Bo;

% Preallocate memory; Initialize time-dependent sum of molecule numbers
time = zeros(1,t_steps);
At = zeros(reps,t_steps);          Bt = zeros(reps,t_steps);          Ct = zeros(reps,t_steps);

extt = zeros(reps, Bo);     % *** To monitor extinction time of B agents ***

for n=1:reps

    fprintf('.');
    
    % Arrays for tracking reactant and product agent status
    A = ones(1,Ao);                 B = ones(1,Bo);                 C = zeros(1,Cmax); 
    At(n,1) = sum(A);               Bt(n,1) = sum(B);               Ct(n,1) = sum(C);

    for t=2:t_steps
        tempB = find(B==1);                 
        P = 1 - exp(-k * At(n,t-1) * dt);

        for i = 1:size(tempB,2)
            if rand < P
                B(tempB(i)) = 0;                    
                tempA = find(A==1);                 x = ceil(rand*size(tempA,2));
                A(tempA(x)) = 0;     
                C(tempB(i)) = 1;    
                
                extt(n,tempB(i)) = time(t-1) + dt; 
            end
        end 

        time(t) = time(t-1) + dt;
        At(n,t) = sum(A);             Bt(n,t) = sum(B);              Ct(n,t) = sum(C);
    end

end                 % end "for reps" loop

%% Population statistics of time trajectories
avgA = mean(At);            avgB = mean(Bt);            avgC = mean(Ct);
sdevA = std(At);            sdevB = std(Bt);            sdevC = std(Ct);

%% Solve ODE for 2nd order kinetics
tmax=time(end);
[t_sol,y_sol] = ode45(@o2_dif,0:tmax/100:tmax,[Ao ; Bo ; Co]);

%% Plot Time Trajectories
figure('Name','2nd Order Rx Time course','NumberTitle','off'); 
plot(time,avgA,time,avgB,time,avgC);                                    hold on;
p1_devA1 = plot(time,avgA+sdevA,'LineStyle','--','Color',[0.8 0.8 0.8]);
p1_devA0 = plot(time,avgA-sdevA,'LineStyle','--','Color',[0.8 0.8 0.8]);
p1_devB1 = plot(time,avgB+sdevB,'LineStyle','--','Color',[0.8 0.8 0.8]);
p1_devB0 = plot(time,avgB-sdevB,'LineStyle','--','Color',[0.8 0.8 0.8]);
p1_devC1 = plot(time,avgC+sdevC,'LineStyle','--','Color',[0.8 0.8 0.8]);
p1_devC0 = plot(time,avgC-sdevC,'LineStyle','--','Color',[0.8 0.8 0.8]);

xlabel('time');             legend('Agent A','Agent B','Agent C');
axis([0 maxTime 0 agents]);

% Plot DE Time Courses
plot(t_sol,y_sol(:,1),':b');               % scatter(ty,y_sol(:,1),'.b');                  
plot(t_sol,y_sol(:,2),':g');               % scatter(ty,y_sol(:,2),'.g');
plot(t_sol,y_sol(:,3),':r');               % scatter(ty,y_sol(:,3),'.r');

clear p1*;

%% Compare simulation standard deviation with CME...
% Obtained from the Chemical Master Equation (see notes in root folder).
% Var = < n^2 > - < n >^2  = < n >           (True in this case; see derivation)
% SDev = sqrt(Var) = sqrt(<n>) 

% Load data from file 'CME_A=Ao_B=Bo.mat'
% Make sure the CME calculations have already been done.     ***********
if exist('sigma','var')==0
    cd ../CME_data;
    load(['CME_A=' num2str(Ao) '_B=' num2str(Bo) '.mat'],'sigma');                       
    cd ../P_ber;
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
SSR_s = sum((sdevB - sigma).^2);                % sum of square residuals (sim data vs DE predictions)

Rsq_s = 1 - SSR_s ./ SST_s                      % Definition of R^2

%% Remove zeros from extt (zeros represent agents that didn't "die" within totalTime)
f = size(extt,2);                      
for z=1:f
    extt_zeros = find(extt(:,z)==0);
    for a = 1:size(extt_zeros,1)
        extt(extt_zeros(a),z) = NaN;
    end
end

extt = reshape(extt,[],1);            % Reshape to column vector

%% Plot Extinction time distribution
bins = 20;

fig1 = figure('WindowStyle','Normal','Position',[1 1 500 450]);
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.

% hist(extt(:,:),bins);                         % plot histogram

[num1, c1] = hist(extt,bins);                   % histogram data 
Num1 = num1 / (reps*Bo);                        % Normalize

p_b = bar(c1',Num1','DisplayName','$\textrm{ABK}$');                               hold on;
set(p_b(1),'FaceColor','r','EdgeColor',[0.70 0.78 1]);        

timeD = 0:maxTime/bins:maxTime;
[t_sol1,y_sol1] = ode45(@o2_dif,timeD,[Ao ; Bo ; Co]);

% dA = abs(diff(y_sol1(:,1))) / Ao;                     
dB = abs(diff(y_sol1(:,2))) / Bo;

timeDm = maxTime/bins/2:maxTime/bins:maxTime;

p_de = plot(timeDm,dB','LineStyle','none','DisplayName','$\textrm{DE}$',...
    'Marker','o','MarkerEdgeColor',[0 0.5 0],'MarkerFaceColor',[0 0.5 0],...
    'MarkerSize',2);

title(['$N_{A,i} = ' num2str(Ao) '\, , \, N_{B,i} = ' num2str(Bo) '$'],...
    'Interpreter','LateX','FontSize',12);

set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('$t \, \textrm{(sec)}$','Interpreter','LateX','FontSize',11);          
ylabel('$\textrm{Fraction of } A+B \rightarrow X \textrm{ Transitions}$',...
    'Interpreter','LateX','FontSize',11);

leg1 = legend([p_b p_de]);
set(leg1,'Interpreter','LateX','FontSize',10,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');       

% Create textbox
annotation(gcf,'textbox',...
    [0.00771 0.93471 0.05010 0.06429],'String','b)',...
    'FontSize',12,'FontName','Times New Roman',...
    'FitBoxToText','off','EdgeColor','none');

% Coefficient of Determination R^2 for histogram fit to distribution PDF.
disp(['R^2 = ' num2str(CoefDet(Num1,dB'))]);

%% Finish
clear temp* a b i f n x w z P;
toc
