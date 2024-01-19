%          R <---->  Rp             Rates: k_f, k_r (Forward rx activated by Yp)
%               ^     |             Michaelis-Menten constants: Km_f, Km_r
%              /      |
%             |       |
%    Y <---> Yp       |             Rates: k_f1, k_r1 (Forward rx activated by X)
%         ^           |             Michaelis-Menten constants: Km_f1, Km_r1
%        /            |
%        |            \
%   -->  X  ---------------->                      Rp: Response
%   k_b  ^       k_d1/k_d2
%        | 
%        | k_s
%        S                             S: Signal

% Simulating negative feedback process (3-component system: X, Y/Yp, R/Rp) 
% k_b: X synthesis, 0th order process
% k_s: X synthesis, 1st order process wrt S (but S is not consummed)
% k_s: ALTERNATIVE IMPLEMENTATION: X synthesis, 0th order process, UPregulated by S
% k_d1: X degradation, 1st order process wrt X
% k_d2: X degradation, 2nd order process wrt X, Rp (but Rp is not consummed)
% k_f: Rp synthesis, MM process wrt Yp (but Yp is not consummed)
% k_r: R synthesis, MM process 
% Km_f, Km_r: MM constants for forward and reverse rxs R <==> Rp
% k_f1: Yp synthesis, MM process wrt X (but X is not consummed)
% k_r1: Y synthesis, MM process 
% Km_f1, Km_r1: MM constants for forward and reverse rxs Y <==> Yp

% Note: Assume number of S molecules/agents is NOT changing (eg, S is DNA molecules)
% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear;       
clc;     tic;                                          
rng(1);
%% Declare Parameters
global agents k_s k_b k_d1 k_d2 k_f k_r Km_f Km_r k_f1 k_r1 Km_f1 Km_r1 S;

reps = 100;                         % Repeat simulation this many times

totalTime = 5000;                   % Simulation time (sec)
dt = 0.02;                          % Constant (fixed) time step increment (sec)
t_steps = totalTime / dt;
time = 0:dt:totalTime;

agents = 250;
k_b = 0;                            % 0th order X synthesis rate
k_d1 = 0;                           % 1st order degradation rate (units: 1/sec)
k_d2 = 0.001;                       % 2nd order degradation rate wrt X, Rp (units: 1/sec)
k_s = 0.40;                         % 0th order R synthesis rate UPregulated by S

k_f = 0.01;                         % basal forward rate R-->Rp
k_r = 1;                            % basal reverse rate R<--Rp
Km_f = 10;                          % MICROSCOPIC Michaelis-Menten constant for forward rx R-->Rp
Km_r = 10;                          % MICROSCOPIC Michaelis-Menten constant for reverse rx R<--Rp

k_f1 = 0.01;                        % basal forward rate Y-->Yp
k_r1 = 1;                           % basal reverse rate Y<--Yp
Km_f1 = 10;                         % MICROSCOPIC Michaelis-Menten constant for forward rx Y-->Yp
Km_r1 = 10;                         % MICROSCOPIC Michaelis-Menten constant for reverse rx Y<--Yp

S = 50;                             % Assume number of S molecules/agents is NOT changing

% Initial population sizes
Xi = 0;
Ypi = 100;                  Yi = agents - Ypi;
Rpi= 100;                   Ri = agents - Rpi;

% Preallocate memory
X_all = zeros(reps,t_steps+1);         
Yp_all = zeros(reps,t_steps+1);        
Rp_all = zeros(reps,t_steps+1);        

%% Solve DE - Get steady state values (assuming damped oscillations)
maxDEtime = 30000;            % in sec, to ensure that steady-states have been reached.
[t_sol, y_sol] = ode15s(@negFb3_dif2,0:maxDEtime/500:maxDEtime,[Xi; Yi; Ypi; Ri; Rpi]);

X_ss = y_sol(end,1);
Y_ss = y_sol(end,2);               Yp_ss = y_sol(end,3);
R_ss = y_sol(end,4);               Rp_ss = y_sol(end,5);

% Symbolic calculations - Linearize and calculate eigenvalues
syms X Yp Rp positive;
dX_sym = + k_b + k_s * S - k_d1 * X - k_d2 * X * Rp;
dYp_sym = + k_f1 * (agents - Yp)   * X / (Km_f1 + (agents - Yp)) - k_r1 * Yp / (Km_r1 + Yp);
dRp_sym = + k_f  * (agents - Rp) * Yp / (Km_f + (agents - Rp)) - k_r * Rp / (Km_r + Rp);

fprintf('Steady-States: \n');
Jac = jacobian([dX_sym, dYp_sym, dRp_sym],[X, Yp, Rp]);             % Calculate Jacobian matrix
J = double(subs(Jac,[X, Yp, Rp],[X_ss, Yp_ss, Rp_ss]));
lambdas = eig(J);
fprintf(['X=' num2str(X_ss) ', Yp=' num2str(Yp_ss) ', Rp=' num2str(Rp_ss) '\n']);
fprintf('Eigenvalues: \n');
disp(num2str(lambdas));                                   clear temp tempX g;

%% Repeat simulation 'reps' number of times
for n=1:reps
        
    fprintf('.');
    
    % Initialize - Preallocate memory
    Tx = zeros(1,t_steps);
    Ty = zeros(1,t_steps);              Typ = zeros(1,t_steps);
    Tr = zeros(1,t_steps);              Trp = zeros(1,t_steps);
    P_f1 = zeros(1,t_steps);            P_r1 = zeros(1,t_steps);
    P_f = zeros(1,t_steps);             P_r = zeros(1,t_steps);
    P_d2 = zeros(1,t_steps);

    % ******** Initial conditions - Number of R, Ep, E Agents ********
    Tx(1) = Xi;    
    Typ(1) = Ypi;      Ty(1) = agents - Typ(1);					
    Trp(1) = Rpi;      Tr(1) = agents - Trp(1);					
    % ****************************************************************

    tempX = zeros(1,2*agents);        
    tempY = zeros(1,agents);        tempR = zeros(1,agents);             
    % Put "1" where agents are "alive", then randomize the array
    for c=1:Tx(1),                  tempX(c)=1;             end
    for d=1:Ty(1),                  tempY(d)=1;             end
    for e=1:Tr(1),                  tempR(e)=1;             end
    tempX = RandArray(tempX);                   % Randomize X array
    tempY = RandArray(tempY);                   % Randomize Y array
    tempR = RandArray(tempR);                   % Randomize R array
    % Markov process, so only previous and current time steps needed --> 2 rows:    
    Xv = [tempX ; tempX];                        % Initialize Vector storing state of X agents
    Yv = [tempY ; tempY];                        % Initialize Vector storing state of Y agents
    Ypv = ~ Yv;                                  % Initialize Vector storing state of Yp agents
    Rv = [tempR ; tempR];                        % Initialize Vector storing state of R agents
    Rpv = ~ Rv;                                  % Initialize Vector storing state of Rp agents
    % - R and Rp are complementary (R + Rp = agents). So are Y and Yp.
    clear c d e tempX tempY tempR;
        
% ABK simulation
    t = 1;

    P_b = k_b * dt;                  % Probability of X synthesis process (0th order)
    P_s = k_s * S * dt;              % Probability of X synthesis process (0th order), UPreg'd by S
    P_d1 = 1 - exp(-k_d1 * dt);      % Probability of X degradation process (1st order) wrt X

    while t <= t_steps % && Tep(t)<agents
        P_d2(t) = k_d2 * Trp(t) * dt;                          % wrt each X molecule
        P_f1(t) = k_f1 * Tx(t) / (Km_f1 + Ty(t)) * dt;         % wrt each Y molecule
        P_r1(t) = k_r1 / (Km_r1 + Typ(t)) * dt;                % wrt each Yp molecule
        P_f(t)  = k_f * Typ(t) / (Km_f + Tr(t)) * dt;          % wrt each R molecule
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

        tempY = find(Yv(1,:)==1);
        for h=1:size(tempY,2)
            if rand < P_f1(t)
                Yv(2,tempY(h)) = 0;         % Conversion of Y ...
                Ypv(2,tempY(h)) = 1;        % to Yp
            end
        end

        tempYp = find(Ypv(1,:)==1);
        for q=1:size(tempYp,2)
            if rand < P_r1(t)
                Ypv(2,tempYp(q)) = 0;       % Conversion of Yp ...
                Yv(2,tempYp(q)) = 1;        % to Y
            end
        end

        tempR = find(Rv(1,:)==1);
        for j=1:size(tempR,2)
            if rand < P_f(t)
                Rv(2,tempR(j)) = 0;         % Conversion of R ...
                Rpv(2,tempR(j)) = 1;        % to Rp
            end
        end

        tempRp = find(Rpv(1,:)==1);
        for k=1:size(tempRp,2)
            if rand < P_r(t)
                Rpv(2,tempRp(k)) = 0;       % Conversion of Rp ...
                Rv(2,tempRp(k)) = 1;        % to R
            end
        end
    
        Tx(t+1) = sum(Xv(2,:));
        Ty(t+1) = sum(Yv(2,:));             Typ(t+1) = sum(Ypv(2,:));          
        Tr(t+1) = sum(Rv(2,:));             Trp(t+1) = sum(Rpv(2,:));      

        Xv(1,:) = Xv(2,:);
        Yv(1,:) = Yv(2,:);                  Ypv(1,:) = Ypv(2,:);    
        Rv(1,:) = Rv(2,:);                  Rpv(1,:) = Rpv(2,:);
        t = t + 1;
    end

    X_all(n,:) = Tx;                Yp_all(n,:) = Typ;                Rp_all(n,:) = Trp;
        
end             % *** end 'for n=1:reps' loop ***

%% Calculate AVERAGE + SDEV Time Course

avg_X = mean(X_all);                    sdev_X = std(X_all);
avg_Yp = mean(Yp_all);                  sdev_Yp = std(Yp_all);
avg_Rp = mean(Rp_all);                  sdev_Rp = std(Rp_all);

N_max = max([max(avg_X+sdev_X) max(avg_Yp+sdev_Yp) max(avg_Rp+sdev_Rp)]);

%% Plot [selected] time course
% Using Report Nomenclature: X = R , Y = W, R = X  (script name = report name) in legends.

trial = 5;

fig0 = figure('Name',['Trial TC, S=' num2str(S) ': SEPARATE Plots'],...
    'NumberTitle','off','Position',[1 1 1000 1218]);        

subplot(2,2,1);                                                                 hold on;
plot(time(1:20:end),X_all(trial,1:20:end),'b','LineWidth',1); 
plot(time(1:20:end),Yp_all(trial,1:20:end),'g','LineWidth',1);
plot(time(1:20:end),Rp_all(trial,1:20:end),'r','LineWidth',1);
% axis([0 time(end) 0 1.05*N_max]);                                               % axis tight; 
axis([0 time(end) 0 250]);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('t (sec)');            ylabel('N(t)');                                        
leg0 = legend('ABK N_R(t)','ABK N_{Wp}(t)','ABK N_{Xp}(t)','Location','NorthEast');
set(leg0,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95]);                                              hold off;

subplot(2,2,2);                                                                 hold on;
plot(time(1:20:end),X_all(trial,1:20:end),'b','LineWidth',1); 
plot(t_sol,y_sol(:,1),'--','Color',[0.15 1 0.75],'LineWidth',2); 
% axis([0 time(end) 0 1.05*N_max]);                                               % axis tight; 
axis([0 time(end) 0 250]);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('t (sec)');            ylabel('N_R(t)');                                        
leg0b = legend('ABK N_R(t)','DE  N_R(t)','Location','NorthEast');
set(leg0b,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95]);                                              hold off;

subplot(2,2,3);                                                                 hold on;
plot(time(1:20:end),Yp_all(trial,1:20:end),'g','LineWidth',1);
plot(t_sol,y_sol(:,3),'--','Color',[0 0.5 0],'LineWidth',2); 
% axis([0 time(end) 0 1.05*N_max]);                                               % axis tight; 
axis([0 time(end) 0 250]);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('t (sec)');            ylabel('N_{Wp}(t)');                                        
leg0c = legend('ABK N_{Wp}(t)','DE  N_{Wp}(t)','Location','NorthEast');
set(leg0c,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95]);                                              hold off;

subplot(2,2,4);                                                                 hold on;
plot(time(1:20:end),Rp_all(trial,1:20:end),'r','LineWidth',1);
plot(t_sol,y_sol(:,5),'--','Color',[0.75 0 1],'LineWidth',2);
% axis([0 time(end) 0 1.05*N_max]);                                               % axis tight; 
axis([0 time(end) 0 250]);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('t (sec)');            ylabel('N_{Xp}(t)');                                        
leg0d = legend('ABK N_{Xp}(t)','DE  N_{Xp}(t)','Location','NorthEast');
set(leg0d,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95]);                                              hold off;

% Create textboxes
annotation(fig0,'textbox',[0.09 0.900 0.05 0.07],'String',{'a)'},'FontSize',12,...
    'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');
annotation(fig0,'textbox',[0.53 0.900 0.05 0.07],'String',{'b)'},'FontSize',12,...
    'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');
annotation(fig0,'textbox',[0.09 0.425 0.05 0.07],'String',{'c)'},'FontSize',12,...
    'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');
annotation(fig0,'textbox',[0.53 0.425 0.05 0.07],'String',{'d)'},'FontSize',12,...
    'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');

%% Plot AVERAGE Time Course - SINGLE plot
% Using Report Nomenclature: X = R , Y = W, R = X  (script name = report name) in legends.

fig1a = figure('Name',['Avg TC, S=' num2str(S) ': SINGLE Plot'],...
    'NumberTitle','off','Position',[1 1 500 406]);                        hold on;

plot(time(1:20:end),avg_X(1:20:end),'b','LineWidth',2);                                          
plot(time(1:20:end),avg_Yp(1:20:end),'g','LineWidth',2);                                
plot(time(1:20:end),avg_Rp(1:20:end),'r','LineWidth',2);

plot(t_sol,y_sol(:,1),'--','Color',[0.15 1 0.75],'LineWidth',2);   
plot(t_sol,y_sol(:,3),'--','Color',[0 0.5 0],'LineWidth',2); 
plot(t_sol,y_sol(:,5),'--','Color',[0.75 0 1],'LineWidth',2);  

axis([0 time(end) 60 250]);                                               
% axis tight;    
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');                   hold off;
xlabel('t (sec)');            ylabel('N(t)');                                        

%% Plot AVERAGE Time Course - SEPARATE plots
% Using Report Nomenclature: X = R , Y = W, R = X  (script name = report name) in legends.

fig1 = figure('Name',['Avg TC, S=' num2str(S) ': SEPARATE Plots'],...
    'NumberTitle','off','Position',[1 1 1000 1218]);        

subplot(2,2,1);                                                                 hold on;
plot(time(1:20:end),avg_X(1:20:end),'b','LineWidth',2);                                          
plot(time(1:20:end),avg_Yp(1:20:end),'g','LineWidth',2);                                
plot(time(1:20:end),avg_Rp(1:20:end),'r','LineWidth',2);
% axis([0 time(end) 0 1.05*N_max]);                                               % axis tight; 
axis([0 time(end) 0 250]);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('t (sec)');            ylabel('N(t)');                                        
leg1 = legend('<N_R(t)>_{sim}','<N_{Wp}(t)>_{sim}','<N_{Xp}(t)>_{sim}','Location','NorthEast');
set(leg1,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95]);                                              hold off;

subplot(2,2,2);                                                                 hold on;
plot(time(1:20:end),avg_X(1:20:end),'b','LineWidth',2);
plot(t_sol,y_sol(:,1),'--','Color',[0.15 1 0.75],'LineWidth',2);   
plot(time(1:20:end),avg_X(1:20:end)+sdev_X(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
plot(time(1:20:end),avg_X(1:20:end)-sdev_X(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
% axis([0 time(end) 0 1.05*N_max]);                                               % axis tight;  
axis([0 time(end) 0 250]);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('t (sec)');            ylabel('N_R(t)');                                        
leg2 = legend('<N_R(t)>_{sim}','DE  N_R(t)','Location','NorthEast');
set(leg2,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95]);                                              hold off;

subplot(2,2,3);                                                                 hold on;
plot(time(1:20:end),avg_Yp(1:20:end),'g','LineWidth',2);         
plot(t_sol,y_sol(:,3),'--','Color',[0 0.5 0],'LineWidth',2); 
plot(time(1:20:end),avg_Yp(1:20:end)+sdev_Yp(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
plot(time(1:20:end),avg_Yp(1:20:end)-sdev_Yp(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
% axis([0 time(end) 0 1.05*N_max]);                                               % axis tight; 
axis([0 time(end) 0 250]);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('t (sec)');            ylabel('N_{Wp}(t)');                                        
leg3 = legend('<N_{Wp}(t)>_{sim}','DE  N_{Wp}(t)','Location','NorthEast');
set(leg3,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95]);                                              hold off;

subplot(2,2,4);                                                                 hold on;
plot(time(1:20:end),avg_Rp(1:20:end),'r','LineWidth',2);
plot(t_sol,y_sol(:,5),'--','Color',[0.75 0 1],'LineWidth',2);  
plot(time(1:20:end),avg_Rp(1:20:end)+sdev_Rp(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
plot(time(1:20:end),avg_Rp(1:20:end)-sdev_Rp(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
% axis([0 time(end) 0 1.05*N_max]);                                             % axis tight; 
axis([0 time(end) 0 250]);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('t (sec)');            ylabel('N_{Xp}(t)');                                        
leg4 = legend('<N_{Xp}(t)>_{sim}','DE  N_{Xp}(t)','Location','NorthEast');
set(leg4,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95]);                                              hold off;

% Create textboxes
annotation(fig1,'textbox',[0.09 0.900 0.05 0.07],'String',{'a)'},'FontSize',12,...
    'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');
annotation(fig1,'textbox',[0.53 0.900 0.05 0.07],'String',{'b)'},'FontSize',12,...
    'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');
annotation(fig1,'textbox',[0.09 0.425 0.05 0.07],'String',{'c)'},'FontSize',12,...
    'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');
annotation(fig1,'textbox',[0.53 0.425 0.05 0.07],'String',{'d)'},'FontSize',12,...
    'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');

%% ****** Spectral Analysis ******
% fs = 1/dt;                                          % sampling frequency
% t_osc = 200:dt:3000;                                % choose time interval with oscillations               
% X_osc = avg_X(t_osc(1)/dt:t_osc(end)/dt);
% Yp_osc = avg_Yp(t_osc(1)/dt:t_osc(end)/dt);                
% Rp_osc = avg_Rp(t_osc(1)/dt:t_osc(end)/dt);
% 
% figure;     plot(t_osc,X_osc,t_osc,Yp_osc,t_osc,Rp_osc);
% xlabel = ('time (sec)');        ylabel('#');
% 
% Osc = [X_osc; Yp_osc; Rp_osc];
f_th = imag(lambdas(2)) / (2*pi);                   % Theoretical oscillation frequency (Hz)
period_th = 1 / f_th;                               % Theoretical period (sec)
% 
% figure('Position',[1 1 350 690]);                   % Make new figure       
% % *** Do FFT on X, Yp, Rp data ***
% n = 2^nextpow2(length(X_osc));                      
% f = fs / n * (0:n/2);                   species = {'X', 'Yp', 'Rp'}; 
% for b=1:3
%     Z(b,:) = fft(Osc(b,:),n);
%     Zn(b,:) = abs(Z(b,:)/n);
%     subplot(3,1,b);                     title(species{b});                      hold on;
%     plot(f,Zn(b,1:n/2+1)); 
%     tempy = [0 25 50];                      tempf = [f_th f_th f_th];           
%     plot(tempf,tempy,'--','Color',[1 0.75 0.25]); 
%     axis([0.5*f_th 1.5*f_th 0 40]);
%     xlabel('f (Hz)');     
%     ylabel('Power');                    hold off;
% end
% clear b tempy tempf species;
%% Finding the period of the first few observed oscillations (avg)
% Adjust time entries accordingly

% Min1 = find(avg_X==min(avg_X(400/dt:600/dt)))*dt;
% Min2 = find(avg_X==min(avg_X(900/dt:1200/dt)))*dt;
% period = Min2(2) - Min1;
% Min3 = find(avg_X==min(avg_X(1400/dt:1700/dt)))*dt;
% period(2) = Min3(5) - Min2(2);
% Min4 = find(avg_X==min(avg_X(1900/dt:2200/dt)))*dt;
% period(3) = Min4(7) - Min3(5)
%% Finish
fprintf(['\n' ElapsedTime(toc)]);

clear array X Rp dX dRp X_nc Rp_nc X_nc_sym Rp_nc_sym Yp leg* fig* T* P* Yv Ypv Rv Rpv Xv;
clear m h i j k n q r t J Jac lambdas temp tempR tempRp tempX tempY tempYp *_sym ans;  
