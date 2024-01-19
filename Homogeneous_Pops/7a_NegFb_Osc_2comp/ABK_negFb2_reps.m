%      R <---->  Rp             Rates: k_f, k_r (Forward rx activated by Yp)
%           ^     |             Michaelis-Menten constants: Km_f, Km_r
%          /      |
%         /       |
%        /        |
%        |        \
%   -->  X  ----------->                      Rp: Response
%   k_b  ^    k_d1/k_d2
%        | 
%        | k_s
%        S                             S: Signal

% Simulating negative feedback process (2-component system: X, R/Rp) 
% k_b: X synthesis, 0th order process
% k_s: X synthesis, 1st order process wrt S (but S is not consummed)
% k_s: ALTERNATIVE IMPLEMENTATION: X synthesis, 0th order process, UPregulated by S
% k_d1: X degradation, 1st order process wrt X
% k_d2: X degradation, 2nd order process wrt X, Rp (but Rp is not consummed)
% k_f: Rp synthesis, MM process wrt Yp (but Yp is not consummed)
% k_r: R synthesis, MM process 
% Km_f, Km_r: MM constants for forward and reverse rxs R <==> Rp

% Note: Assume number of S molecules/agents is NOT changing (eg, S is DNA molecules)
% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear;       
clc;     tic;                                          
rng(1);
%% Declare Parameters
global agents k_s k_b k_d1 k_d2 k_f k_r Km_f Km_r S;

reps = 100;                         % Repeat simulation this many times

totalTime = 2500;                   % Simulation time (sec)
dt = 0.01;                          % Constant (fixed) time step increment (sec)
t_steps = totalTime / dt;
time = 0:dt:totalTime;

agents = 1000;
k_b = 0.00;                         % 0th order X synthesis rate
k_d1 = 0.00;                        % 1st order degradation rate (units: 1/sec)
k_d2 = 0.001;                       % 2nd order degradation rate wrt X, Rp (units: 1/sec)
k_s = 0.30;                         % 0th order R synthesis rate UPregulated by S
k_f = 0.0015;                       % basal forward rate (1st order)
k_r = 0.8;                          % basal reverse rate (1st order)
Km_f = 30;                          % MICROSCOPIC Michaelis-Menten constant for forward rx
Km_r = 1;                           % MICROSCOPIC Michaelis-Menten constant for reverse rx

S = 30;                             % Assume number of S molecules/agents is NOT changing

% Initial population sizes
Xi = 0;                Rpi= 0;            Ri = agents - Rpi;

% Preallocate memory
X_all = zeros(reps,t_steps+1);         
Rp_all = zeros(reps,t_steps+1);        

%% Solve DE - Get steady state values (assuming damped oscillations)
maxDEtime = 2500;            % in sec, to ensure that steady-states have been reached.
[t_sol, y_sol] = ode45(@negFb_dif,0:maxDEtime/500:maxDEtime,[Xi; Rpi]);

X_ss = y_sol(end,1);                Rp_ss = y_sol(end,2);

% Symbolic calculations - Linearize and calculate eigenvalues
syms X Rp positive;
dX_sym = + k_b + k_s * S - k_d1 * X - k_d2 * X * Rp;
dRp_sym = + k_f  * (agents - Rp) * X / (Km_f + (agents-Rp)) - k_r * Rp / (Km_r + Rp);

fprintf('Steady-States: \n');
Jac = jacobian([dX_sym,dRp_sym],[X, Rp]);             % Calculate Jacobian matrix
J = double(subs(Jac,[X, Rp],[X_ss, Rp_ss]));
lambdas = eig(J);
fprintf(['X=' num2str(X_ss) ', Rp=' num2str(Rp_ss) '\n']);
fprintf('Eigenvalues: \n');
disp(num2str(lambdas));                                   clear temp tempX g;

%% Repeat simulation 
for n=1:reps
        
    fprintf('.');
    
    % Initialize, set Initial conditions.
    Tx = zeros(1,t_steps);
    Tr = zeros(1,t_steps);              Trp = zeros(1,t_steps);
    P_f = zeros(1,t_steps);             P_r = zeros(1,t_steps);
    P_d2 = zeros(1,t_steps);

    % ******** Initial conditions - Number of R, Ep, E Agents ********
    Tx(1) = Xi;    
    Trp(1) = Rpi;                       Tr(1) = agents - Trp(1);					
    % ****************************************************************

    tempX = zeros(1,1.5*agents);          tempR = zeros(1,agents);             
    % Put "1" where agents are "alive", then randomize the array
    for c=1:Tx(1),                      tempX(c)=1;             end
    for e=1:Tr(1),                      tempR(e)=1;             end
    tempX = RandArray(tempX);                    % Randomize X array
    tempR = RandArray(tempR);                    % Randomize R array
%   Markov process, so only previous and current time steps needed --> 2 rows:    
    Xv = [tempX ; tempX];                        % Initialize Vector storing state of X agents
    Rv = [tempR ; tempR];                        % Initialize Vector storing state of R agents
    Rpv = ~ Rv;                                  % Initialize Vector storing state of Rp agents
    % - R and Rp are complementary (R + Rp = agents).
    clear c d e tempX tempR;
        
% ABK simulation
    t = 1;

    P_b = k_b * dt;                  % Probability of X synthesis process (0th order)
    P_s = k_s * S * dt;              % Probability of X synthesis process (0th order), UPreg'd by S
    P_d1 = 1 - exp(-k_d1 * dt);      % Probability of X degradation process (1st order) wrt X

    while t <= t_steps % && Tep(t)<agents
        P_d2(t) = k_d2 * Trp(t) * dt;                          % wrt each X molecule
        P_f(t)  = k_f * Tx(t) / (Km_f + Tr(t)) * dt;          % wrt each R molecule
        P_r(t)  = k_r / (Km_r + Trp(t)) * dt;                  % wrt each Rp molecule

        % Take care of 0th order processes first
        if rand < P_b
            tempX = find(Xv(1,:)==0);                % Randomly choose X agent synthesis 
            Xv(2,tempX(ceil(rand * size(tempX,2)))) = 1;
        end

        if rand < P_s
            tempX = find(Xv(1,:)==0);                % Randomly choose X agent synthesis 
            Xv(2,tempX(ceil(rand * size(tempX,2)))) = 1;
        end
        % End of 0th order processes

    % The following implementation also works (treating it as a 1st order rx wrt S)     
    %     for h=1:S                               % Reaction for each S molecule/agent
    %         if rand < P_s
    %             tempX = find(Xv(1,:)==0);       % Randomly choose X agent which becomes "alive" 
    %             Xv(2,tempR(ceil(rand * size(tempX,2)))) = 1;
    %         end
    %     end

        P_d = P_d1 + P_d2(t);               % Total probability of X degradation
        tempX = find(Xv(1,:)==1);
        for i = 1:size(tempX,2)
            if rand < P_d                   % Degradation of X, 2 contributing processes
                Xv(2,tempX(i)) = 0;         % X agent is degraded
            end
        end

        tempR = find(Rv(1,:)==1);
        for j=1:size(tempR,2)
            if rand < P_f(t)
                Rv(2,tempR(j)) = 0;         % Conversion of R to Rp
                Rpv(2,tempR(j)) = 1;
            end
        end

        tempRp = find(Rpv(1,:)==1);
        for k=1:size(tempRp,2)
            if rand < P_r(t)
                Rpv(2,tempRp(k)) = 0;       % Conversion of Rp to R
                Rv(2,tempRp(k)) = 1;
            end
        end

        Tx(t+1) = sum(Xv(2,:));
        Tr(t+1) = sum(Rv(2,:));             Trp(t+1) = sum(Rpv(2,:));      

        Xv(1,:) = Xv(2,:);
        Rv(1,:) = Rv(2,:);                  Rpv(1,:) = Rpv(2,:);
        t = t + 1;
    end

    X_all(n,:) = Tx;                      Rp_all(n,:) = Trp;
    
%     clear Tx Tr Trp Xv Rv Rpv;
%     clear P_b P_s P_d P_d1 P_d2 P_f P_r;
    
end             % *** end 'for n=1:reps' loop ***

%% Calculate AVERAGE , SDEV Time Course, and Coefficient of Variation

avg_X = mean(X_all);                 sdev_X = std(X_all);           cv_X = sdev_X ./ avg_X;       
avg_Rp = mean(Rp_all);               sdev_Rp = std(Rp_all);         cv_Rp = sdev_Rp ./ avg_Rp;  

N_max = max([max(avg_X+sdev_X) max(avg_Rp+sdev_Rp)]);

%% Plot [selected] time course
% Using Report Nomenclature: X -> R , R -> X  (script name -> report name) in legend.

trial = 5;                         % Choose simulation run/trial (trial <= reps)

fig0 = figure('Name','Time Course','NumberTitle','off');                    hold on;
set(fig0,'Position',[1 1 500 450]);
fig0.PaperUnits = 'inches';
fig0.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.

plot(time(1:20:end),X_all(trial,1:20:end),'b'); 
plot(time(1:20:end),Rp_all(trial,1:20:end),'r');
plot(t_sol,y_sol(:,1),'--','Color',[0.15 1 0.75],'LineWidth',2); 
plot(t_sol,y_sol(:,2),'--','Color',[0.75 0 1],'LineWidth',2);  

axis([0 time(end) 0 N_max+130]);                                               hold off;
xlabel('t (sec)');                 ylabel('N(t)');   

% leg0 = legend('ABK N_R(t)','ABK N_{Xp}(t)','DE N_R(t)','DE N_{Xp}(t)','Location','Best');
leg0 = legend('ABK N_R(t)','ABK N_{Xp}(t)','Location','Best');            % Without DEs
set(leg0,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');

% savefig('TC_sample.fig');
clear trial;
%% Plot AVERAGE Time Course
% Using Report Nomenclature: X = R , R = X  (script name = report name) in legend.

fig1 = figure('Name','Avg TC','NumberTitle','off');                              hold on;
set(fig1,'Position',[1 1 500 450]);
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.

plot(time(1:20:end),avg_X(1:20:end),'b','LineWidth',2);
plot(time(1:20:end),avg_Rp(1:20:end),'r','LineWidth',2);
plot(t_sol,y_sol(:,1),'--','Color',[0.15 1 0.75],'LineWidth',2); 
plot(t_sol,y_sol(:,2),'--','Color',[0.75 0 1],'LineWidth',2);  
plot(time(1:20:end),avg_X(1:20:end)+sdev_X(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
plot(time(1:20:end),avg_X(1:20:end)-sdev_X(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
plot(time(1:20:end),avg_Rp(1:20:end)+sdev_Rp(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
plot(time(1:20:end),avg_Rp(1:20:end)-sdev_Rp(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);

axis([0 time(end) 0 1.05*N_max]);                                               % axis tight;    
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('t (sec)');            ylabel('N(t)');                                        
leg1 = legend('ABK <N_R(t)>','ABK <N_{Xp}(t)>','DE  N_R(t)','DE  N_{Xp}(t)');
set(leg1,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');                       hold off;

% savefig('TC_avg.fig');
%% Plot coefficient of variation
fig2 = figure('Name','Coefficient of Variation','NumberTitle','off');                  
set(fig2,'Position',[1 1 500 450]);                                     hold on;                  
fig2.PaperUnits = 'inches';
fig2.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.

p1 = plot(time,cv_X,'k','DisplayName','ABK, R','LineWidth',2);  
% pt1 = plot(t_sol,1./sqrt(y_sol(:,1)),':','Color',[0.2 0.75 0.5],...
%     'DisplayName','Poisson, R','LineWidth',2);
pt1 = plot(time,1./sqrt(avg_X),':','Color',[0.2 0.75 0.5],...
    'DisplayName','Poisson, R','LineWidth',2);

p2 = plot(time,cv_Rp,'Color',[0.5 0.5 0.5],'DisplayName','ABK, Xp','LineWidth',2);               
% pt2 = plot(t_sol,1./sqrt(y_sol(:,2)),':','Color',[0.1 0.5 0.75],...
%     'DisplayName','Poisson, Xp','LineWidth',2);
pt2 = plot(time,1./sqrt(avg_Rp),':','Color',[0.1 0.5 0.75],...
    'DisplayName','Poisson, Xp','LineWidth',2);

% set([p(2) p(3) pt(2) pt(3)],'Visible','off');       % Don't show for ...
axis([0 totalTime 0 1]);        
% axis tight;
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');

leg2 = legend([p1 pt1 p2 pt2]);
set(leg2,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95],'Location','Best');

xlabel('t (sec)');         
ylabel('$\eta$','Interpreter','LaTeX');                                 hold off;                                    

% savefig('CV.fig');
%% ****** Spectral Analysis ******
% fs = 1/dt;                                          % sampling frequency
% t_osc = 200:dt:3000;                                % choose time interval with oscillations               
% X_osc = avg_X(t_osc(1)/dt:t_osc(end)/dt);
% Rp_osc = avg_Rp(t_osc(1)/dt:t_osc(end)/dt);
% 
% figure;     plot(t_osc,X_osc,t_osc,Rp_osc);
% xlabel = ('time (sec)');        ylabel('#');
% 
% Osc = [X_osc; Rp_osc];
% f_th = imag(lambdas(2)) / (2*pi);                   % Theoretical oscillation frequency (Hz)
% period_th = 1 / f_th;                               % Theoretical period (sec)
% 
% figure('Position',[1 1 350 690]);                   % Make new figure       
% % *** Do FFT on X, Yp, Rp data ***
% n = 2^nextpow2(length(X_osc));                      
% f = fs / n * (0:n/2);                   species = {'X', 'Yp', 'Rp'}; 
% for b=1:2
%     Z(b,:) = fft(Osc(b,:),n);
%     Zn(b,:) = abs(Z(b,:)/n);
%     subplot(3,1,b);                     title(species{b});                      hold on;
%     plot(f,Zn(b,1:n/2+1)); 
%     tempy = [0 25 50];                      tempf = [f_th f_th f_th];           
%     plot(tempf,tempy,'--','Color',[1 0.75 0.25]); 
%     axis([0.5*f_th 1.5*f_th 0 40]);
%     xlabel('f (Hz)');     % gives ERROR: Index exceeds matrix dimensions. Weird. 
%     ylabel('Power');                    hold off;
% end
% clear b tempy tempf species;

%% Finish
fprintf(['\n' ElapsedTime(toc)]);

clear array X Rp dX dRp X_nc Rp_nc X_nc_sym Rp_nc_sym Yp leg* fig* T* P* Yv Ypv Rv Rpv Xv;
clear m h i j k n p* q r t J Jac lambdas temp tempR tempRp tempX tempY tempYp *_sym ans;  
