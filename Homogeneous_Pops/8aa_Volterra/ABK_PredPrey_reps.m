% Volterra Predator-Prey system
%     R --> 2R              Density-dependent Birth of R; rate constant: a 
% F + R -->  F              Death of R; rate constant: b (2nd order)
% F     -->                 Death of F; rate constant: c (1st order)
% F + R --> 2F + R          Birth of F; rate constant: d (2nd order)

% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

% Note two (limiting) cases that come up:
% - If R goes extinct, then F will too (according to this model).
% - If F goes extinct, then R grows exponentially.

clear;      clc;     tic;                           
rng(1);

global a b c d K;

reps = 500;                    % Repeat simulation this many times
% BUT, I will discard simulations where the population of R blows up
% (when F goes extinct), or both species go extinct.
to_remove = [];

totalTime = 100;                 % Simulation time (sec)
dt = 0.01;                      % Constant (fixed) time step increment (sec)
t_steps = totalTime / dt + 1;
time = 0:dt:totalTime;

agents = 300;                  % Max number of agents allowed
% If a population size exceeds this, then species F is extinct and R grows exponentially.

a = 1;                          % a rate constant (1st order) 
b = 0.01;                       % b rate constant (2nd order)
c = 1;                          % c rate constant (1st order) 
d = 0.01;                       % d rate constant (2nd order)
K = 500;                        % Carrying capacity for R

% *** Initial population sizes ***
Ri = 100;                        Fi = 100;

% *** Set up differential equations of system using symbolic variables ***
syms R F positive;
dR_sym = a*R*(1-R/K) - b*F*R;
dF_sym = d*R*F - c*F;

R_ss = [0 c/d]                     
F_ss = [0 a/b * (1 - c/(d*K))]

Jac = jacobian([dR_sym,dF_sym],[R,F]);

% Linearize and calculate eigenvalues
for w=1:size(R_ss,2)   
    J = double(subs(Jac,[R,F],[R_ss(w),F_ss(w)]));
    TrJ = J(1,1) + J(2,2);
    DetJ = det(J);
    DiscrJ = TrJ^2 - 4*DetJ;
    lambda_plus(w) = (TrJ + sqrt(DiscrJ)) / 2;
    lambda_minus(w) = (TrJ - sqrt(DiscrJ)) / 2;        
end
clear w;

disp([num2str(lambda_plus)]);               disp([num2str(lambda_minus)]);

% Conditions for having a Stable Spiral
cond1 = 4*d*K * (d*K/c - 1);
if a < cond1,    disp('Stable Spiral predicted (wrt param a).');   end

cond2 = ( 2*d/a * (sqrt(1 + a/c) - 1) )^-1;
if K > cond2,    disp('Stable Spiral predicted (wrt param K).');   end

% ************** End symbolic calculations (for now) *****************

% * Preallocate memory *
R_all = zeros(reps,t_steps);          F_all = zeros(reps,t_steps);         

%% Repeat simulation 'reps' number of times
for n=1:reps
        
    fprintf('.');
    
    % Initialize - Preallocate memory
    Tr = zeros(1,t_steps);              Tf = zeros(1,t_steps);

    Pb = zeros(1,t_steps);              Pd = zeros(1,t_steps);

    % ******** Initial conditions - Number of R, F Agents ************
    Tr(1) = Ri;                         Tf(1) = Fi;    
    % ****************************************************************

    tempR = zeros(1,agents);            tempF = zeros(1,agents);                   
    % Put "1" where agents are "alive", then randomize the array
    for x=1:Tr(1),                      tempR(x)=1;             end
    tempR = RandArray(tempR);           % Randomize array
    Rv = [tempR ; tempR];               % Initialize vector storing state of R agents     

    for y=1:Tf(1),                      tempF(y)=1;             end
    tempF = RandArray(tempF);           % Randomize array
    Fv = [tempF ; tempF];               % Initialize vector storing state of F agents 
    % Note on Rv, Fv: Markov process, so only previous and current time steps needed --> 2 rows:          
    clear x y tempR tempF;

    % *** ABK simulation ***
    t = 1;

    Pc = 1 - exp(-c * dt);                      % Pr for death of F: F --> 

    while t < t_steps && Tr(t) <= K && Tf(t) ~= 0

        Pa(t) = a * (1 - Tr(t)/K) * dt;         % Pr for birth of R: R --> 2R
        Pb(t) = b * Tf(t) * dt;                 % Pr for death of R: F + R --> F  (wrt R)
        Pd(t) = d * Tr(t) * dt;                 % Pr for birth of F: F + R --> 2F + R  (wrt F)

        tempR_al = find(Rv(1,:)==1);            tempR_dd = find(Rv(1,:)==0);
        
        for i = 1:size(tempR_al,2)
            r = rand;
            if r < Pb(t)
                Rv(2,tempR_al(i)) = 0;         % R(i) agent dies     (2nd order)         
            elseif r >= Pb(t) && r < Pb(t)+Pa(t)   % R agent 'i' gives birth (1st order)
                j = ceil(rand * size(tempR_dd,2));   % to randomly chosen R agent 'j'
                Rv(2,tempR_dd(j)) = 1;         % R(j) agent is born
            end
        end

        tempF_al = find(Fv(1,:)==1);            tempF_dd = find(Fv(1,:)==0);
        
        for k = 1:size(tempF_al,2)
            r = rand;
            if r < Pc                               % **check probability condition**
                Fv(2,tempF_al(k)) = 0;                 % F(k) agent dies
            elseif r >= Pc && r < Pc+Pd(t)       % F agent 'm' is born (2nd order)
                m = ceil(rand * size(tempF_dd,2));     % randomly chose F agent 'm'
                Fv(2,tempF_dd(m)) = 1;                 % F(m) agent is born
            end
        end
        
        Tr(t+1) = sum(Rv(2,:));              Tf(t+1) = sum(Fv(2,:));          

        Rv(1,:) = Rv(2,:);                   Fv(1,:) = Fv(2,:);

        t = t + 1;
    end         % end "while" loop

    if t < t_steps
        to_remove = [to_remove n];
    end
    
    R_all(n,:) = Tr;                F_all(n,:) = Tf;   

end         % end "for n=1:reps" loop

clear Rv Fv;
%%  Remove runs with extinctions (both R, F) or R blowup (F extinction)
for z = size(to_remove,2):-1:1
    R_all(to_remove(z),:) = [];
    F_all(to_remove(z),:) = [];
end

extinction_free_reps = size(R_all,1)
%% Calculate AVERAGE + SDEV Time Course (START HERE WHEN LOADING DATA FILE)
avg_R = mean(R_all);                    sdev_R = std(R_all);            cv_R = sdev_R ./ avg_R; 
avg_F = mean(F_all);                    sdev_F = std(F_all);            cv_F = sdev_F ./ avg_F; 
%% Solve DE
% if exist('t','var')==0,         finaltime = totalTime;      
% else                            finaltime = (t-1) * dt;         end
[t_sol, y_sol] = ode45(@predprey_dif,0:totalTime/1000:totalTime,[Ri ; Fi]);

%% Plot [selected] time course
trial = 6;
ul = max(R_all(trial,:));

fig0 = figure('Name','Time Course','NumberTitle','off','Position',[1 1 500 450]);      hold on;
fig0.PaperUnits = 'inches';
fig0.PaperPosition = [0 0 6 5];                 % Control the size of printed (eps) fig.

% title(['$N_{R,i} = ' num2str(Ri) '\, , \, N_{F,i} = ' num2str(Fi) '$'],...
%     'FontSize',12,'Interpreter','LateX');

p_Rs = plot(time(1:20:end),R_all(trial,1:20:end),'r','DisplayName','$\textrm{ABK } N_{R}(t)$');                                
p_Fs = plot(time(1:20:end),F_all(trial,1:20:end),'g','DisplayName','$\textrm{ABK } N_{F}(t)$');

p_Rd = plot(t_sol,y_sol(:,1),'--','Color',[0.75 0 1],'LineWidth',2,...
    'DisplayName','$\textrm{DE } N_{R}(t)$');                                       
p_Fd = plot(t_sol,y_sol(:,2),'--','Color',[0 0.5 0],'LineWidth',2,...
    'DisplayName','$\textrm{DE } N_{F}(t)$');

% axis tight;                                     
axis([0 time(end) 0 210]);                               
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');                           hold off;
xlabel('$t \; \textrm{(months)}$','Interpreter','LateX','FontSize',11);                 
ylabel('$N(t)$','Interpreter','LateX','FontSize',11); 

leg0 = legend([p_Rs , p_Fs , p_Rd , p_Fd]);           % ABK + DE trajectories
% leg0 = legend([p_Rs , p_Fs]);                       % Only ABK trajectories
set(leg0,'Interpreter','LateX','FontSize',10,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');       

annotation(gcf,'textbox',...
    [0.00771 0.93471 0.05010 0.06429],'String','a)',...
    'FontSize',12,'FontName','Times New Roman',...
    'FitBoxToText','off','EdgeColor','none');

%% Spectral Analysis of 'trial' trajectory for both species
fs = 1/dt;                                          % sampling frequency
t_osc_start = 40;
% t_osc = t_osc_start:dt:totalTime; 
R_osc = R_all(trial,t_osc_start/dt+1:end);
F_osc = F_all(trial,t_osc_start/dt+1:end);          
% plot(t_osc,R_osc,t_osc,F_osc);     xlabel = ('time (sec)');        ylabel('#');

Osc = [R_osc; F_osc];
f_th = abs(imag(lambda_plus(2))) / (2*pi);          % Theoretical oscillation frequency (Hz)
period_th = 1 / f_th;                               % Theoretical period (sec)

fig1 = figure('Name',['Spectral Analysis: trial #' num2str(trial)],...
    'NumberTitle','off','Position',[501 1 350 550]);
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0 0 6 10];                 % Control the size of printed (eps) fig.

% do FFT on R, F data
n = 2^(nextpow2(length(F_osc))+4);                      
f = fs / n * (0:n/2);                   

species = {'\underbar{Species R}', '\underbar{Species F}'}; 

Z = zeros(size(species,2),n);               Zn = zeros(size(species,2),n);

for q=1:size(species,2)
    Z(q,:) = fft(Osc(q,:),n);
    Zn(q,:) = abs(Z(q,:)/n);
    subplot(size(species,2),1,q);                              hold on;
    title(species{q},'FontName','Times New Roman','FontSize',12,'Interpreter','LaTex');  
    plot(1000*f,Zn(q,1:n/2+1));                                   
    
%     axis([0.5*f_th 2*f_th 0 12*q]);
    xlim([100*f_th 1800*f_th]);
    ylim([0 2*Zn(q,find(abs(f-f_th)<=0.002 & abs(f-f_th)>=0.0005,1))]);

    tempy = [0 50 100];                     tempf = 1000*[f_th f_th f_th];           
    plot(tempf,tempy,'r');                  
    set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
    xlabel('f (mHz)');                      ylabel('Power');
end

% Create textbox
annotation(fig1,'textbox',[0.03 0.92 0.05 0.07],'String',{'b)'},'FontSize',12,...
    'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');

% hold off;                               
% clear q tempy tempf species;

%% Plot AVERAGE Time Course - SINGLE PLOT

fig2 = figure('Name','Avg TC','NumberTitle','off','Position',[801 1 500 450]);      hold on;
fig2.PaperUnits = 'inches';
fig2.PaperPosition = [0 0 6 5];                 % Control the size of printed (eps) fig.

% title(['$N_{R,i} = ' num2str(Ri) '\, , \, N_{F,i} = ' num2str(Fi) '$'],...
%     'FontSize',12,'Interpreter','LateX');                            

p_R = plot(time(1:20:end),avg_R(1:20:end),'r','LineWidth',2,...
    'DisplayName','$\textrm{ABK } < \! N_{R}(t) \! >$');   
plot(time(1:20:end),avg_R(1:20:end)+sdev_R(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
plot(time(1:20:end),avg_R(1:20:end)-sdev_R(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);

p_F = plot(time(1:20:end),avg_F(1:20:end),'g','LineWidth',2,...
    'DisplayName','$\textrm{ABK } < \! N_{F}(t) \! >$');
% plot(time(1:20:end),avg_F(1:20:end)+sdev_F(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
% plot(time(1:20:end),avg_F(1:20:end)-sdev_F(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);

p_Rd = plot(t_sol,y_sol(:,1),'--','Color',[0.75 0 1],'LineWidth',2,...
    'DisplayName','$\textrm{DE } N_{R}(t)$');
p_Fd = plot(t_sol,y_sol(:,2),'--','Color',[0 0.5 0],'LineWidth',2,...
    'DisplayName','$\textrm{DE } N_{F}(t)$');                                     
% alpha(p_Rd,0.5);
% alpha(p_Fd,0.5);

% axis tight;                                         
axis([0 time(end) 0 210]);                                                   
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');                   hold off;
xlabel('$t \; \textrm{(months)}$','Interpreter','LateX','FontSize',11);            
ylabel('$N(t)$','Interpreter','LateX','FontSize',11);                                        

leg2 = legend([p_R , p_F , p_Rd , p_Fd]);
set(leg2,'Interpreter','LateX','FontSize',10,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');                                          

annotation(gcf,'textbox',...
    [0.00771 0.93471 0.05010 0.06429],'String','b)',...
    'FontSize',12,'FontName','Times New Roman',...
    'FitBoxToText','off','EdgeColor','none');

%% Plot AVERAGE Time Course - SEPARATE plots
fig1 = figure('Name','Avg TC ; SEPARATE Plots',...
    'NumberTitle','off','Position',[1 1 1000 500]);        

subplot(1,2,1);                                                        hold on;         
plot(time(1:20:end),avg_R(1:20:end),'r','LineWidth',3);
plot(t_sol,y_sol(:,1),'--','Color',[0.75 0 1],'LineWidth',2);   
plot(time(1:20:end),avg_R(1:20:end)+sdev_R(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
plot(time(1:20:end),avg_R(1:20:end)-sdev_R(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
% axis([0 time(end) 0 1.05*N_max]);                                 
% axis tight;  
axis([0 time(end) 0 200]);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('t (sec)');            ylabel('N_R(t)');                                        
leg2 = legend('ABK <N_R(t)>','DE  N_R(t)','Location','NorthEast');
set(leg2,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95]);                                     hold off;

subplot(1,2,2);                                                        hold on;
plot(time(1:20:end),avg_F(1:20:end),'g','LineWidth',3);
plot(t_sol,y_sol(:,2),'--','Color',[0 0.5 0],'LineWidth',2);  
plot(time(1:20:end),avg_F(1:20:end)+sdev_F(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
plot(time(1:20:end),avg_F(1:20:end)-sdev_F(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
% axis([0 time(end) 0 1.05*N_max]);                                % axis tight; 
axis([0 time(end) 0 200]);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('t (sec)');            ylabel('N_{X}(t)');                                        
leg3 = legend('ABK <N_{X}(t)>','DE  N_{X}(t)','Location','NorthEast');
set(leg3,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95]);                                     hold off;

% Create textboxes
annotation(fig1,'textbox',[0.09 0.900 0.05 0.07],'String',{'a)'},'FontSize',12,...
    'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');
annotation(fig1,'textbox',[0.09 0.425 0.05 0.07],'String',{'b)'},'FontSize',12,...
    'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');

%% Spectral Analysis of the AVG trajectories for both species
fs = 1/dt;                                          % sampling frequency
t_osc = 0:dt:totalTime;                          
R_osc = avg_R(t_osc(1)/dt+1:end);
F_osc = avg_F(t_osc(1)/dt+1:end);

% plot(t_osc,R_osc,t_osc,F_osc);      % xlabel = ('time (sec)');      ylabel('#');

Osc = [R_osc; F_osc];
% f_th = abs(imag(lambda_plus(2))) / (2*pi);          % Theoretical oscillation frequency (Hz)
% period_th = 1 / f_th;                               % Theoretical period (sec)

fig3 = figure('Name','Spectral Analysis: AVG','NumberTitle','off','Position',[801 1 350 550]);  
fig3.PaperUnits = 'inches';
fig3.PaperPosition = [0 0 6 5];                 % Control the size of printed (eps) fig.

% do FFT on R, Ep, X data
n = 2^(nextpow2(length(F_osc))+4);                      
f = fs / n * (0:n/2);                   species = {'Species R', 'Species F'};  

Z  = zeros(size(species,2),n);
Zn = zeros(size(species,2),n);

for q=1:size(species,2)
    Z(q,:) = fft(Osc(q,:),n);
    Zn(q,:) = abs(Z(q,:)/n);
    subplot(size(species,2),1,q);        title(species{q});              hold on;
    plot(1000*f,Zn(q,1:n/2+1));  
    
%     axis([0.5*f_th 2*f_th 0 12*q]);
    xlim([100*f_th 1800*f_th]);
    ylim([0 6*Zn(q,find(abs(f-f_th)<=0.002 & abs(f-f_th)>=0.0005,1))]);


    tempy = [0 50 100];                     tempf = 1000*[f_th f_th f_th];           
    plot(tempf,tempy,'r');                  
    set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
    xlabel('f (mHz)');                      ylabel('Power');
end

hold off;                                   clear q tempy tempf species;

%% Plot Nullclines, State-Space trajectories
% ** Nullclines ** (explicitly stated)
R_nc = a/b .* (1 - [0:5:agents]/K);     % a/K * R^2 + R*(bF-a) = 0  implies F = a/b * (1 - R/K)               
F_nc = c/d;                             % R = c/d;
 
[Rm,Fm] = meshgrid(40:8:160);     % Mesh-grid values for constructing state space
dR_num = double(subs(dR_sym,[R,F],{Rm,Fm}));
dF_num = double(subs(dF_sym,[R,F],{Rm,Fm}));
% r = ( dR_num.^2 + dF_num.^2 ).^0.5;
r=1;

fig4 = figure('Name','State Space','NumberTitle','off','Position',[1 1 500 450]);          hold on;   
fig4.PaperUnits = 'inches';
fig4.PaperPosition = [0 0 6 5];                 % Control the size of printed (eps) fig.

p_ncR = plot(0:5:agents,R_nc,'Color',[1 0.3 0.2],...
    'LineWidth',2,'DisplayName','$N_R \textrm{ nullcline}$');
p_ncF = plot([F_nc F_nc],[0 160],'Color',[0.4 0.4 1],...
    'LineWidth',2,'DisplayName','$N_F \textrm{ nullcline}$');            

% p_df = quiver(Rm,Fm,dR_num./r,dF_num./r,'g');                    % Direction field  
p_df = streakarrow(Rm,Fm,dR_num./r,dF_num./r,0.7,1);             % Direction field

p_de = plot(y_sol(:,1),y_sol(:,2),'k','MarkerSize',2,...
    'DisplayName','$\textrm{DE trajectory}$');                 % plot deterministic trajectory
% p_abk = plot(avg_R,avg_F,'b','LineWidth',1,...
%     'DisplayName','$\textrm{Mean ABK trajectory}$');            % plot stochastic trajectory

p_ic = plot(Ri,Fi,'rp','MarkerSize',9,'MarkerFaceColor','r',...
    'DisplayName',['$\textrm{Initial: } (' num2str(Ri) '\,,' num2str(Fi) ')$']);     % plot initial condition

% StableFP_name = ['$\textrm{Stable FP (Spiral): } (' num2str(R_ss(2)) '\,,' num2str(F_ss(2)) ')$'];
StableFP_name = '$\textrm{Stable FP (Spiral)}$';

p_fp = plot(R_ss(2),F_ss(2),'oc','MarkerSize',8,...
    'MarkerFaceColor','b',...
    'DisplayName',StableFP_name);  % FP

set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
axis([60 140 40 120]);          set(gca,'DataAspectRatio',[1 1 1]);            
xlabel('$N_R$','Interpreter','LateX','FontSize',11);                  
ylabel('$N_F$','Interpreter','LateX','FontSize',11);                  

leg4 = legend([p_ncR, p_ncF, p_fp, p_ic, p_de]);
set(leg4,'Interpreter','LateX','FontSize',10,...   
    'EdgeColor',[0.95 0.95 0.95]);                                                  hold off;

%% AVG Power spectrum of ALL trajectories for species R: Calculation
% Warning: Memory intensive calculation ( > 4GB)
fs = 1/dt;                                          % sampling frequency
f_th = abs(imag(lambda_plus(2))) / (2*pi);          % Theoretical oscillation frequency (Hz)
period_th = 1 / f_th;                               % Theoretical period (sec)

% do FFT on R data
n = 2^(nextpow2(size(R_all,2))+4);                      
f = fs / n * (0:n/2);                   

t_osc_start = 25;

Z = zeros(reps,n);               Zn = zeros(reps,n);

for q=1:reps
    Z(q,:) = fft(R_all(q,t_osc_start/dt+1:end),n);
    Zn(q,:) = abs(Z(q,:)/n);
end

Zn_avg = mean(Zn);                             % clear Z Zn q;
%% Now plot avg power spectrum
fig5 = figure('Name','AVG Power spectrum of all trials',...
    'NumberTitle','off','Position',[501 1 500 450]);
fig5.PaperUnits = 'inches';                                     hold on;
fig5.PaperPosition = [0 0 6 5];                 % Control the size of printed (eps) fig.

% title(species{q},'FontName','Times New Roman','FontSize',12,'Interpreter','LaTeX');  
plot(1000*f,Zn_avg(1:n/2+1),'LineWidth',2);                                   

tempy = [0 50 100];                     tempf = 1000*[f_th f_th f_th];    
plot(tempf,tempy,'r');         

%     axis([0.5*f_th 2*f_th 0 12*q]);
xlim([350*f_th 1500*f_th]);
ylim([0 1.2*Zn(q,find(abs(f-f_th)<=0.002 & abs(f-f_th)>=0.0005,1))]);

set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('$f \times 10^{-3} \, \textrm{month}^{-1}$','Interpreter','LateX','FontSize',11);                      
ylabel('$\textrm{Power}$','Interpreter','LateX','FontSize',11);                     hold off;

annotation(gcf,'textbox',...
    [0.00771 0.93471 0.05010 0.06429],'String','d)',...
    'FontSize',12,'FontName','Times New Roman',...
    'FitBoxToText','off','EdgeColor','none');

clear tempy tempf;
%% Plot Coefficient of variation
fig6 = figure('Name','CV',...
    'NumberTitle','off','Position',[501 1 500 450]);
fig6.PaperUnits = 'inches';                                     hold on;
fig6.PaperPosition = [0 0 6 5];                 % Control the size of printed (eps) fig.

p1 = plot(time,cv_R,'k','DisplayName','$\textrm{ABK, } R$','LineWidth',2);  
% pt1 = plot(t_sol,1./sqrt(y_sol(:,1)),':','Color',[0.2 0.75 0.5],...
%     'DisplayName','Poisson, R','LineWidth',2);
pt1 = plot(time,1./sqrt(avg_R),':','Color',[0.2 0.75 0.5],...
    'DisplayName','$\textrm{Poisson, } R$','LineWidth',2);

p2 = plot(time,cv_F,'Color',[0.5 0.5 0.5],'DisplayName','$\textrm{ABK, } F$','LineWidth',2);               
% pt2 = plot(t_sol,1./sqrt(y_sol(:,2)),':','Color',[0.1 0.5 0.75],...
%     'DisplayName','Poisson, Xp','LineWidth',2);
pt2 = plot(time,1./sqrt(avg_F),':','Color',[0.1 0.5 0.75],...
    'DisplayName','$\textrm{Poisson, } F$','LineWidth',2);

% set([p(2) p(3) pt(2) pt(3)],'Visible','off');       % Don't show for ...
axis([0 totalTime 0 0.6]);                                        hold off;  
% axis tight;
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');

leg6 = legend([p1 pt1 p2 pt2]);
set(leg6,'Interpreter','LateX','FontSize',10,...
    'EdgeColor',[0.95 0.95 0.95],'Location','Best');

xlabel('$t \; \textrm{(months)}$','Interpreter','LateX','FontSize',11);         
ylabel('$\eta$','Interpreter','LaTeX','FontSize',11);                                                                     

annotation(gcf,'textbox',...
    [0.00771 0.93471 0.05010 0.06429],'String','c)',...
    'FontSize',12,'FontName','Times New Roman',...
    'FitBoxToText','off','EdgeColor','none');

%% Confirm if E is a conserved quantity for deterministic solution
% ****** Conserved quantity E(R,F) ?? ********

% E_sym = c*log(R) + a*log(F) - d*R - b*F;
% 
% for v=1:size(y_sol,1)
%     E(v) = double(subs(E_sym,[R F],[y_sol(v,1) y_sol(v,2)]));
% end
% plot(diff(E));

% Result: E is NOT conserved! (as expected for the Volterra model)
%% Finish
fprintf(['\n' ElapsedTime(toc)]);

clear R Ep dR* dF* R_nc F_nc;
clear h i j k m n r t v temp* tempR* texpF*;
clear fig* leg* p_*;
