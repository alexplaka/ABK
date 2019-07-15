%          Ep  <--->  E          Rates: k_f, k_r (Reverse rx activated by R)
%          |     ^               Michaelis-Menten constants: Km_f, Km_r
%           \     \              k_s: synthesis of X
%    k_s     \     |             k_b1,2: synthesis of R
%  S ---> X  --->  R  ------>    k_d: degradation of R                
%           k_b1,2      k_d               ** S: Signal, R: Response **  

% Simulating substrate depletion oscillator process: 
% k_s: X synthesis, 0th order process, UPregulated by S
% k_b1: X --> R, 1st order process wrt X
% k_b2: X --> R, 2nd order process wrt X, Ep (Ep is not consummed)
% k_d: R degradation, 1st order process wrt R
% k_f: E synthesis, MM process
% k_r: Ep synthesis, MM process wrt R (but R is not consummed)
% Km_f, Km_r: MM constants for forward and reverse rxs

% Note: Assume number of S molecules/agents is NOT changing (eg, S is DNA molecules)

% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;      clc;     tic;                           % rng(1);

global agents k_s k_b1 k_b2 k_d k_f k_r Km_f Km_r S;

reps = 500;                         % Repeat simulation this many times

totalTime = 5000;                   % Simulation time (sec)
dt = 0.02;                          % Constant (fixed) time step increment (sec)
t_steps = totalTime / dt;
time = 0:dt:totalTime;

agents = 100;
k_s = 0.075;                         % X synthesis, 0th order process, UPregulated by S
k_b1 = 0.001;                       % X --> R, 1st order process wrt X
k_b2 = 0.003;                       % X --> R, 2nd order process wrt X, Ep (Ep is not consummed)
k_d = 0.04;                         % R degradation, 1st order process wrt R
k_f = 0.5;                        % Ep --> E forward rate
k_r = 0.0125;                        % Ep <-- E reverse rate, UPregulated by R
Km_f = 5;                      % MICROSCOPIC Michaelis-Menten constant for forward rx
Km_r = 5;                      % MICROSCOPIC Michaelis-Menten constant for reverse rx
S = 10;                            % Assume number of S molecules/agents is NOT changing

% Initial population sizes
Ri = 0;                 Epi = 0;                Xi = 0;

% Preallocate memory
X_all = zeros(reps,t_steps+1);         
Ep_all = zeros(reps,t_steps+1);        
R_all = zeros(reps,t_steps+1);        

%% Repeat simulation 'reps' number of times
for n=1:reps
        
    fprintf('.');
    
    % Initialize - Preallocate memory
    Tr = zeros(1,t_steps);              Tx = zeros(1,t_steps);
    Tep = zeros(1,t_steps);             Te = zeros(1,t_steps);

    P_f = zeros(1,t_steps);             P_r = zeros(1,t_steps);
    P_b2 = zeros(1,t_steps);            

    % ******** Initial conditions - Number of R, Ep, E Agents ********
    Tr(1) = Ri;                       Tx(1) = Xi;    
    Tep(1) = Epi;                     Te(1) = agents - Epi;					
    % ****************************************************************

    tempR = zeros(1,2*agents);        tempEp = zeros(1,agents); 
    tempX = zeros(1,3*agents);
    % Put "1" where agents are "alive", then randomize the array
    for c=1:Tr(1),                  tempR(c)=1;             end
    for d=1:Tep(1),                 tempEp(d)=1;            end
    for e=1:Tx(1),                  tempX(e)=1;             end
    tempR = RandArray(tempR);                   % Randomize R array
    tempEp = RandArray(tempEp);                 % Randomize Ep array
    tempX = RandArray(tempX);                   % Randomize X array
    % Markov process, so only previous and current time steps needed --> 2 rows:    
    R = [tempR ; tempR];                        % Initialize vector storing state of R agents
    Ep = [tempEp ; tempEp];                     % Initialize vector storing state of Ep agents
    E = ~ Ep;                                   % Initialize vector storing state of E agents
    % - Ep and E are complementary (E + Ep = agents)
    X = [tempX ; tempX];                        % Initialize vector storing state of X agents
    clear c d e tempR tempEp tempX;

    % *** ABK simulation ***
    t = 1;

    P_s = k_s * S * dt;             % Probability of X synthesis process (0th order)
    P_b1 = k_b1 * dt;               % Probability X-->R (1st order wrt X)
    P_d = k_d * dt;                 % Probability of R degradation process (1st order wrt R)

    while t <= t_steps % && Tep(t)<agents

        P_b2(t) = k_b2 * Tep(t) * dt;   % Probability (wrt X) of R synthesis (2nd order)
        P_f(t) = k_f / (Km_f + Tep(t)) * dt;          % wrt each Ep molecule
        P_r(t) = k_r * Tr(t) / (Km_r + Te(t)) * dt;   % wrt each E molecule

        % Take care of 0th order processes first    
        if rand < P_s
            tempX = find(X(1,:)==0);    % Randomly choose X agent synthesis 
            X(2,tempX(ceil(rand * size(tempX,2)))) = 1;
        end
        % End of 0th order processes

        tempR = find(R(1,:)==1);
        for i = 1:size(tempR,2)
            if rand < P_d                   % Degradation of R
                R(2,tempR(i)) = 0;          % R agent is degraded
            end
        end

        tempE = find(E(1,:)==1);
        for j=1:size(tempE,2)
            if rand < P_r(t)
                E(2,tempE(j)) = 0;          % Conversion of E to Ep
                tempEp = find(Ep(1,:)==0);  % Randomly choose Ep agent synthesis
                Ep(2,tempEp(ceil(rand * size(tempEp,2)))) = 1;
            end
        end

        tempEp = find(Ep(1,:)==1);
        for k=1:size(tempEp,2)
            if rand < P_f(t)
                Ep(2,tempEp(k)) = 0;        % Conversion of Ep to E
                tempE = find(E(1,:)==0);    % Randomly choose E agent synthesis
                E(2,tempE(ceil(rand * size(tempE,2)))) = 1;
            end
        end

        P_b = P_b1 + P_b2(t);               % Total Probability of X --> R
        tempX = find(X(1,:)==1);
        for m = 1:size(tempX,2)
            if rand < P_b                   % X-->R, from 2 separate processes
                X(2,tempX(m)) = 0;          % X agent is converted to R
                tempR = find(R(1,:)==0);    % Randomly choose R agent synthesis 
                R(2,tempR(ceil(rand * size(tempR,2)))) = 1;            
            end
        end  

        Tr(t+1) = sum(R(2,:));              Tx(t+1) = sum(X(2,:));          
        Tep(t+1) = sum(Ep(2,:));            Te(t+1) = sum(E(2,:));

        R(1,:) = R(2,:);                    X(1,:) = X(2,:);
        Ep(1,:) = Ep(2,:);                  E(1,:) = E(2,:);

        t = t + 1;
    end         % end "while" loop

    % Remove unnecessary terminal 0's from arrays
    if t < t_steps
        Tr = Tr(1:t);                               Tx = Tx(1:t);                             
        Tep = Tep(1:t);                             Te = Te(1:t);             
        P_f = P_f(1:t);          P_r = P_r(1:t);    P_b2 = P_b2(1:t);
    end
    
    R_all(n,:) = Tr;                Ep_all(n,:) = Tep;                X_all(n,:) = Tx;

end         % end "for n=1:reps" loop

clear R Ep E X;
%% Calculate AVERAGE + SDEV Time Course

avg_R = mean(R_all);                    sdev_R = std(R_all);            cv_R = sdev_R ./ avg_R;
avg_Ep = mean(Ep_all);                  sdev_Ep = std(Ep_all);          cv_Ep = sdev_Ep ./ avg_Ep;
avg_X = mean(X_all);                    sdev_X = std(X_all);            cv_X = sdev_X ./ avg_X; 

%% Solve DE
if exist('t','var')==0,         finaltime = totalTime;      
else                            finaltime = (t-1) * dt;         end
[t_sol, y_sol] = ode45(@subdep2_dif,0:finaltime/1000:finaltime,[Ri ; Epi ; Xi]);
%% Plot [selected] time course
trial = 4;

figure('Name','Time Course','NumberTitle','off','Position',[1 1 650 406]);           hold on;

p_Rs = plot(time(1:20:end),R_all(trial,1:20:end),'r','DisplayName','ABK N_R(t)');                                
p_Eps = plot(time(1:20:end),Ep_all(trial,1:20:end),'g','DisplayName','ABK N_{Ep}(t)');
% p_Xs = plot(time(1:20:end),X_all(trial,1:20:end),'b','DisplayName','ABK N_X(t)');

p_Rd = plot(t_sol,y_sol(:,1),'--','Color',[0.75 0 1],'LineWidth',2,'DisplayName','DE N_R(t)');                                       
p_Epd = plot(t_sol,y_sol(:,2),'--','Color',[0 0.5 0],'LineWidth',2,'DisplayName','DE N_{Ep}(t)');
% p_Xd = plot(t_sol,y_sol(:,3),'--','Color',[0.15 1 0.75],'LineWidth',2,'DisplayName','DE N_X(t)');                                     

axis tight;                                     % axis([0 finaltime 0 agents]);                               
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');                           hold off;
xlabel('t (sec)');                 ylabel('N(t)'); 

% leg0 = legend([p_Rs , p_Eps , p_Xs , p_Rd , p_Epd , p_Xd]);       % All species
leg0 = legend([p_Rs , p_Rd]);       % Just R
set(leg0,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEastOutside');                                              
%% Symbolic calculations
syms R Ep X positive;
dR_sym = + k_b1 * X + k_b2 * Ep * X - k_d * R;
dEp_sym = + k_r * (agents - Ep) * R / (Km_r + (agents - Ep)) - k_f * Ep / (Km_f + Ep);
dX_sym = + k_s * S - k_b1 * X - k_b2 * Ep * X;

% ** Numerical solution of steady-state values **
SS = vpasolve([dR_sym,dEp_sym,dX_sym],[R,Ep,X]);    % Num solve for equns = 0
R_ss = double(SS.R);
Ep_ss = double(SS.Ep);
X_ss = double(SS.X);
% ***********************************************

fprintf(1,'Steady-States:\t\t\tEigenvalues\n');
Jac = jacobian([dR_sym,dEp_sym,dX_sym],[R,Ep,X]);

for w=1:size(X_ss,1)
    J = double(subs(Jac,[R,Ep,X],[R_ss(w),Ep_ss(w),X_ss(w)]));
    lambdas = eig(J);
    fprintf(['R=' num2str(R_ss(w)) ', Ep=' num2str(Ep_ss(w)) ', X=' num2str(X_ss(w)) ':\n']);
    disp([num2str(lambdas)]);
end
clear w;

%% Spectral Analysis of 'trial' trajectory
fs = 1/dt;                                          % sampling frequency
t_osc = 500:dt:totalTime;                          X_osc = X_all(trial,t_osc(1)/dt+1:end);
Ep_osc = Ep_all(trial,t_osc(1)/dt+1:end);          R_osc = R_all(trial,t_osc(1)/dt+1:end);
% figure;     plot(t_osc,X_osc,t_osc,Yp_osc,t_osc,Rp_osc);
% xlabel = ('time (sec)');        ylabel('#');
Osc = [R_osc; Ep_osc; X_osc];
f_th = abs(imag(lambdas(2))) / (2*pi);              % Theoretical oscillation frequency (Hz)
period_th = 1 / f_th;                               % Theoretical period (sec)

figure('Name',['Spectral Analysis: trial #' num2str(trial)],'NumberTitle','off','Position',[1 1 350 800]);

% do FFT on R, Ep, X data
n = 2^(nextpow2(length(X_osc))+4);                      
f = fs / n * (0:n/2);                   

species = {'\underbar{Species R}', '\underbar{Species Ep}', '\underbar{Species X}'}; 

Z = zeros(size(species,2),n);               Zn = zeros(size(species,2),n);

for b=1:size(species,2)
    Z(b,:) = fft(Osc(b,:),n);
    Zn(b,:) = abs(Z(b,:)/n);
    subplot(3,1,b);                                             hold on;
    title(species{b},'FontName','Times New Roman','FontSize',12,'Interpreter','LaTex');  
    plot(1000*f,Zn(b,1:n/2+1));                                   
    
%     axis([0.5*f_th 2*f_th 0 12*b]);
%     xlim([100*f_th 1800*f_th]);
    xlim([0 6]);                    % Adjust this when Im(lambdas) = 0 
    ylim([0 2*Zn(b,find(abs(f-f_th)<=0.0002 & abs(f-f_th)>=0.00005,1))]);


    tempy = [0 50 100];                     tempf = 1000*[f_th f_th f_th];           
    plot(tempf,tempy,'r');                  
    set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
    xlabel('f (mHz)');                      ylabel('Power');
end
  
hold off;                               clear b tempy tempf species;

%% Plot AVERAGE Time Course - SINGLE PLOT

figure('Name',['Avg TC, S=' num2str(S)],...
    'NumberTitle','off','Position',[1 1 650 406]);                        hold on;

p_R = plot(time(1:20:end),avg_R(1:20:end),'r','LineWidth',2,'DisplayName','<N_R(t)>_{sim}');   
plot(time(1:20:end),avg_R(1:20:end)+sdev_R(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
plot(time(1:20:end),avg_R(1:20:end)-sdev_R(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);

p_Ep = plot(time(1:20:end),avg_Ep(1:20:end),'g','LineWidth',2,'DisplayName','<N_{Ep}(t)>_{sim}'); 
% plot(time(1:20:end),avg_Ep(1:20:end)+sdev_Ep(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
% plot(time(1:20:end),avg_Ep(1:20:end)-sdev_Ep(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);

p_X = plot(time(1:20:end),avg_X(1:20:end),'b','LineWidth',2,'DisplayName','<N_{X}(t)>_{sim}');
% plot(time(1:20:end),avg_X(1:20:end)+sdev_X(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
% plot(time(1:20:end),avg_X(1:20:end)-sdev_X(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);

p_Rd = plot(t_sol,y_sol(:,1),'--','Color',[0.75 0 1],'LineWidth',2,'DisplayName','DE N_R(t)');                                       
p_Epd = plot(t_sol,y_sol(:,2),'--','Color',[0 0.5 0],'LineWidth',2,'DisplayName','DE N_{Ep}(t)');
p_Xd = plot(t_sol,y_sol(:,3),'--','Color',[0.15 1 0.75],'LineWidth',2,'DisplayName','DE N_X(t)');                                     

axis tight;                                         % axis([0 time(end) 0 agents]);                                                   
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');                   hold off;
xlabel('t (sec)');            ylabel('N(t)');                                        

leg1 = legend([p_R , p_Ep , p_X , p_Rd , p_Epd , p_Xd]);
set(leg1,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEastOutside');                                          

%% Plot AVERAGE Time Course - SEPARATE plots

fig1 = figure('Name',['Avg TC, S=' num2str(S) ': SEPARATE Plots'],...
    'NumberTitle','off','Position',[1 1 1000 1218]);        

subplot(2,2,1); 
plot(time(1:20:end),avg_R(1:20:end),'r','LineWidth',2);                 hold on;
plot(time(1:20:end),avg_Ep(1:20:end),'g','LineWidth',2);  
plot(time(1:20:end),avg_X(1:20:end),'b','LineWidth',2);                                          
% axis([0 time(end) 0 1.05*N_max]);                                               
axis tight; 
% axis([0 time(end) 0 250]);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('t (sec)');            ylabel('N(t)');                                        
leg1 = legend('<N_R(t)>_{sim}','<N_{Ep}(t)>_{sim}','<N_{E}(t)>_{sim}','Location','NorthEast');
set(leg1,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95]);                                              hold off;

subplot(2,2,2);                                                                 hold on;
plot(time(1:20:end),avg_R(1:20:end),'r','LineWidth',2);
plot(t_sol,y_sol(:,1),'--','Color',[0.75 0 1],'LineWidth',2);   
plot(time(1:20:end),avg_R(1:20:end)+sdev_R(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
plot(time(1:20:end),avg_R(1:20:end)-sdev_R(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
% axis([0 time(end) 0 1.05*N_max]);                                               % axis tight;  
axis([0 time(end) 0 200]);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('t (sec)');            ylabel('N_R(t)');                                        
leg2 = legend('<N_R(t)>_{sim}','DE  N_R(t)','Location','NorthEast');
set(leg2,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95]);                                              hold off;

subplot(2,2,3);                                                                 hold on;
plot(time(1:20:end),avg_Ep(1:20:end),'g','LineWidth',2);         
plot(t_sol,y_sol(:,2),'--','Color',[0 0.5 0],'LineWidth',2); 
plot(time(1:20:end),avg_Ep(1:20:end)+sdev_Ep(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
plot(time(1:20:end),avg_Ep(1:20:end)-sdev_Ep(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
% axis([0 time(end) 0 1.05*N_max]);                                               % axis tight; 
axis([0 time(end) 0 70]);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('t (sec)');            ylabel('N_{Ep}(t)');                                        
leg3 = legend('<N_{Ep}(t)>_{sim}','DE  N_{Ep}(t)','Location','NorthEast');
set(leg3,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95]);                                              hold off;

subplot(2,2,4);                                                                 hold on;
plot(time(1:20:end),avg_X(1:20:end),'b','LineWidth',2);
plot(t_sol,y_sol(:,3),'--','Color',[0.15 1 0.75],'LineWidth',2);  
plot(time(1:20:end),avg_X(1:20:end)+sdev_X(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
plot(time(1:20:end),avg_X(1:20:end)-sdev_X(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
% axis([0 time(end) 0 1.05*N_max]);                                             % axis tight; 
axis([0 time(end) 0 400]);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('t (sec)');            ylabel('N_{X}(t)');                                        
leg4 = legend('<N_{X}(t)>_{sim}','DE  N_{X}(t)','Location','NorthEast');
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

%% Plot coefficient of variation
fig2 = figure('Name','Coefficient of Variation','NumberTitle','off');                  
fig2.PaperUnits = 'inches';
fig2.PaperPosition = [0 0 6 5];          % Control the size of printed (eps) fig.
set(fig2,'Position',[1 1 500 450]);                                     hold on;                  

p2 = plot(time(1:20:end),cv_Ep(1:20:end),'Color',[0.2 0.75 0.5],...
    'DisplayName','ABK, Ep','LineWidth',2);               
% pt2 = plot(t_sol,1./sqrt(y_sol(:,2)),':','Color',[0.1 0.5 0.75],...
%     'DisplayName','Poisson, Xp','LineWidth',2);
pt2 = plot(time(1:20:end),1./sqrt(avg_Ep(1:20:end)),':','Color',[0.5 0.75 0.2],...
    'DisplayName','Poisson, Ep','LineWidth',2);
p3 = plot(time(1:20:end),cv_X(1:20:end),'Color',[0.35 0.5 0.75],...
    'DisplayName','ABK, X','LineWidth',2);               
pt3 = plot(time(1:100:end),1./sqrt(avg_X(1:100:end)),':','Color',[0 0.5 0.75],...
    'DisplayName','Poisson, X','LineWidth',2);
p1 = plot(time(1:20:end),cv_R(1:20:end),'k','DisplayName','ABK, R','LineWidth',2);               
pt1 = plot(time(1:20:end),1./sqrt(avg_R(1:20:end)),':','Color',[0.5 0.5 0.5],...
    'DisplayName','Poisson, R','LineWidth',2);

% set([p2 p3 pt2 pt3],'Visible','off');           % Don't show for ...
axis([0 totalTime 0 2.5]);        
% axis tight;
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');

leg2 = legend([p1 pt1 p2 pt2 p3 pt3]);              % All in legend
% leg2 = legend([p1 pt1]);                            % for plotting just one species
set(leg2,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');

xlabel('t (sec)');         
ylabel('$\eta$','Interpreter','LaTeX');                                 hold off;                                    

%% State-space 3D 
figure('Name','State Space 3D','NumberTitle','off','Position',[1 1 500 406]);           hold on;

plot3(y_sol(:,3),y_sol(:,2),y_sol(:,1),'k','LineWidth',1.5);        % plot deterministic trajectory
plot3(SS.X(2),SS.Ep(2),SS.R(2),'oc','MarkerSize',8,'MarkerFaceColor','b','DisplayName','Stable FP');                            % plot fixed point

axis([0 150 0 20 0 50]);                                % axis tight;
xlabel('N_X(t)');            ylabel('N_{Ep}(t)');           zlabel('N_R(t)');                   
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');                               hold off;

%% Spectral Analysis of the AVG trajectories
fs = 1/dt;                                          % sampling frequency
t_osc = 500:dt:totalTime;                          X_osc = avg_X(t_osc(1)/dt+1:end);
Ep_osc = avg_Ep(t_osc(1)/dt+1:end);                R_osc = avg_R(t_osc(1)/dt+1:end);

% figure;     plot(t_osc,R_osc,t_osc,Ep_osc,t_osc,X_osc);
% xlabel = ('time (sec)');        ylabel('#');

Osc = [R_osc; Ep_osc; X_osc];
f_th = abs(imag(lambdas(2))) / (2*pi);              % Theoretical oscillation frequency (Hz)
period_th = 1 / f_th;                               % Theoretical period (sec)

figure('Name','Spectral Analysis: AVG','NumberTitle','off','Position',[1 1 350 800]);          

% do FFT on R, Ep, X data
n = 2^(nextpow2(length(X_osc))+4);                      
f = fs / n * (0:n/2);                   species = {'Species R', 'Species Ep', 'Species X'};  

Z  = zeros(size(species,2),n);
Zn = zeros(size(species,2),n);

for b=1:size(species,2)
    Z(b,:) = fft(Osc(b,:),n);
    Zn(b,:) = abs(Z(b,:)/n);
    subplot(3,1,b);        title(species{b});              hold on;
    plot(f,Zn(b,1:n/2+1));  
    
%     axis([0 2*f_th 0 10*b]);
%     xlim([0.1*f_th 2*f_th]);
    xlim([0 0.006]);                    % Adjust this when Im(lambdas) = 0 
    ylim([0 2.5*Zn(b,find(abs(f-f_th)<=0.0002 & abs(f-f_th)>=0.00005,1))]);

    tempy = [0 50 100];                     tempf = [f_th f_th f_th];           
    plot(tempf,tempy,'r');                  
    ylabel('Power');
end
xlabel('f (Hz)');   
hold off;                                   clear b tempy tempf species;

%% Make plots with only species R
plot_R_only;

%% Finish
fprintf(['\n' ElapsedTime(toc)]);

clear array R Ep dR dEp R_nc Ep_nc R_nc_sym Ep_nc_sym;
clear h i j k m n r t temp tempR tempEp tempE tempX;
clear fig0 leg0 fig* leg* p_*;
