% Lotka-Volterra Predator-Prey system 
%     R --> 2R              Birth of R; rate constant: a (1st order)
% F + R -->  F              Death of R; rate constant: b (2nd order)
% F     -->                 Death of F; rate constant: c (1st order)
% F + R --> 2F + R          Birth of F; rate constant: d (2nd order)

% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

% Note two (limiting) cases that come up:
% - If R goes extinct, then F will too (according to this model).
% - If F goes extinct, then R grows exponentially.

clear;      clc;     tic;                           % rng(1);

global a b c d;

reps = 1000;                    % Repeat simulation this many times
% BUT, I will discard simulations where the population of R blows up
% (when F goes extinct), or both species go extinct.
to_remove = [];

totalTime = 50;                 % Simulation time (sec)
dt = 0.01;                      % Constant (fixed) time step increment (sec)
t_steps = totalTime / dt + 1;
time = 0:dt:totalTime;

agents = 1000;                  % Max number of agents allowed
% If a population size exceeds this, then species F is extinct and R grows exponentially.

a = 1;                          % a rate constant (1st order) 
b = 0.01;                       % b rate constant (2nd order)
c = 1;                          % c rate constant (1st order) 
d = 0.01;                       % d rate constant (2nd order)

% *** Initial population sizes ***
Ri = 100;                        Fi = 100;

% *** Set up differential equations of system using symbolic variables ***
syms R F positive;
dR_sym = a*R - b*F*R;
dF_sym = d*R*F - c*F;

R_ss = [0 c/d];                         F_ss = [0 a/b];

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

% ****** Conserved quantity E(R,F) ********
E_sym = c*log(R) + a*log(F) - d*R - b*F;

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

    Pa = 1 - exp(-a * dt);                      % Pr for birth of R: R --> 2R
    Pc = 1 - exp(-c * dt);                      % Pr for death of F: F --> 

    while t < t_steps && Tr(t) < 0.9*agents && Tf(t)~=0

        Pb(t) = b * Tf(t) * dt;                 % Pr for death of R: F + R --> F  (wrt R)
        Pd(t) = d * Tr(t) * dt;                 % Pr for birth of F: F + R --> 2F + R  (wrt F)

        tempR_al = find(Rv(1,:)==1);            tempR_dd = find(Rv(1,:)==0);
        for i = 1:size(tempR_al,2)
            r = rand;
            if r < Pb(t)
                Rv(2,tempR_al(i)) = 0;         % R(i) agent dies     (2nd order)         
            elseif r >= Pb(t) && r < Pb(t)+Pa   % R agent 'i' gives birth (1st order)
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
%% Calculate AVERAGE + SDEV Time Course
avg_R = mean(R_all);                    sdev_R = std(R_all);
avg_F = mean(F_all);                    sdev_F = std(F_all);

%% Solve DE
% if exist('t','var')==0,         finaltime = totalTime;      
% else                            finaltime = (t-1) * dt;         end
[t_sol, y_sol] = ode45(@predprey_dif,0:totalTime/1000:totalTime,[Ri ; Fi]);

%% Plot [selected] time course
trial = 1;
ul = max(R_all(trial,:));

fig0 = figure('Name','Time Course','NumberTitle','off','Position',[1 1 500 500]);      hold on;
title(['N_{R,i} = ' num2str(Ri) ' , N_{F,i} = ' num2str(Fi)],...
    'FontSize',12,'FontName','Times New Roman');

p_Rs = plot(time(1:20:end),R_all(trial,1:20:end),'r','DisplayName','ABK N_{R}(t)');                                
p_Fs = plot(time(1:20:end),F_all(trial,1:20:end),'g','DisplayName','ABK N_{F}(t)');

p_Rd = plot(t_sol,y_sol(:,1),'--','Color',[0.75 0 1],'LineWidth',2,'DisplayName','DE N_{R}(t)');                                       
p_Fd = plot(t_sol,y_sol(:,2),'--','Color',[0 0.5 0],'LineWidth',2,'DisplayName','DE N_{F}(t)');

% axis tight;                                     
axis([0 time(end) 0 ul-2]);                               
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');                           hold off;
xlabel('t (months)');                 ylabel('N(t)'); 

leg0 = legend([p_Rs , p_Fs , p_Rd , p_Fd]);           % ABK + DE trajectories
% leg0 = legend([p_Rs , p_Fs]);                       % Only ABK trajectories
set(leg0,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');       

% Create textbox
annotation(fig0,'textbox',[0.03 0.92 0.05 0.07],'String',{'a)'},'FontSize',12,...
    'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');

%% Spectral Analysis of 'trial' trajectory
fs = 1/dt;                                          % sampling frequency
t_osc = 0:dt:totalTime; 
R_osc = R_all(trial,1:end);
F_osc = F_all(trial,1:end);          
% plot(t_osc,R_osc,t_osc,F_osc);     xlabel = ('time (sec)');        ylabel('#');

Osc = [R_osc; F_osc];
f_th = abs(imag(lambda_plus(2))) / (2*pi);          % Theoretical oscillation frequency (Hz)
period_th = 1 / f_th;                               % Theoretical period (sec)

fig1 = figure('Name',['Spectral Analysis: trial #' num2str(trial)],...
    'NumberTitle','off','Position',[1 1 350 550]);

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
    ylim([0 25*Zn(q,find(abs(f-f_th)<=0.002 & abs(f-f_th)>=0.0005,1))]);


    tempy = [0 50 100];                     tempf = 1000*[f_th f_th f_th];           
    plot(tempf,tempy,'r');                  
    set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
    xlabel('f (mHz)');                      ylabel('Power');
end

% Create textbox
% annotation(fig1,'textbox',[0.03 0.92 0.05 0.07],'String',{'b)'},'FontSize',12,...
%     'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');

% hold off;                               
% clear q tempy tempf species;

%% Plot AVERAGE Time Course - SINGLE PLOT

fig2 = figure('Name','Avg TC','NumberTitle','off','Position',[1 1 500 500]);  
title(['N_{R,i} = ' num2str(Ri) ' , N_{F,i} = ' num2str(Fi)],...
    'FontSize',12,'FontName','Times New Roman');                            hold on;

p_R = plot(time(1:20:end),avg_R(1:20:end),'r','LineWidth',2,'DisplayName','<N_{R}(t)>_{sim}');   
plot(time(1:20:end),avg_R(1:20:end)+sdev_R(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
plot(time(1:20:end),avg_R(1:20:end)-sdev_R(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);

p_F = plot(time(1:20:end),avg_F(1:20:end),'g','LineWidth',2,'DisplayName','<N_{F}(t)>_{sim}');
% plot(time(1:20:end),avg_F(1:20:end)+sdev_F(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
% plot(time(1:20:end),avg_F(1:20:end)-sdev_F(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);

p_Rd = plot(t_sol,y_sol(:,1),'--','Color',[0.75 0 1],'LineWidth',1,'DisplayName','DE N_R(t)');
p_Fd = plot(t_sol,y_sol(:,2),'--','Color',[0 0.5 0],'LineWidth',1,'DisplayName','DE N_{F}(t)');                                     
% alpha(p_Rd,0.5);
% alpha(p_Fd,0.5);

% axis tight;                                         
axis([0 time(end) 0 ul-2]);                                                   
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');                   hold off;
xlabel('t (months)');            ylabel('N(t)');                                        

leg2 = legend([p_R , p_F , p_Rd , p_Fd]);
set(leg2,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');                                          

% Create textbox
annotation(fig2,'textbox',[0.03 0.92 0.05 0.07],'String',{'b)'},'FontSize',12,...
    'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');

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
axis([0 time(end) 0 500]);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('t (sec)');            ylabel('N_R(t)');                                        
leg2 = legend('<N_R(t)>_{sim}','DE  N_R(t)','Location','NorthEast');
set(leg2,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95]);                                     hold off;

subplot(1,2,2);                                                        hold on;
plot(time(1:20:end),avg_F(1:20:end),'g','LineWidth',3);
plot(t_sol,y_sol(:,2),'--','Color',[0 0.5 0],'LineWidth',2);  
plot(time(1:20:end),avg_F(1:20:end)+sdev_F(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
plot(time(1:20:end),avg_F(1:20:end)-sdev_F(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
% axis([0 time(end) 0 1.05*N_max]);                                % axis tight; 
axis([0 time(end) 0 500]);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('t (sec)');            ylabel('N_{X}(t)');                                        
leg3 = legend('<N_{X}(t)>_{sim}','DE  N_{X}(t)','Location','NorthEast');
set(leg3,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95]);                                     hold off;

% Create textboxes
annotation(fig1,'textbox',[0.09 0.900 0.05 0.07],'String',{'a)'},'FontSize',12,...
    'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');
annotation(fig1,'textbox',[0.09 0.425 0.05 0.07],'String',{'b)'},'FontSize',12,...
    'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');

%% Spectral Analysis of the AVG trajectories
fs = 1/dt;                                          % sampling frequency
t_osc = 0:dt:totalTime;                          
R_osc = avg_R(t_osc(1)/dt+1:end);
F_osc = avg_F(t_osc(1)/dt+1:end);

% plot(t_osc,R_osc,t_osc,F_osc);      % xlabel = ('time (sec)');      ylabel('#');

Osc = [R_osc; F_osc];
% f_th = abs(imag(lambda_plus(2))) / (2*pi);          % Theoretical oscillation frequency (Hz)
% period_th = 1 / f_th;                               % Theoretical period (sec)

figure('Name','Spectral Analysis: AVG','NumberTitle','off','Position',[1 1 350 550]);          

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
    ylim([0 25*Zn(q,find(abs(f-f_th)<=0.002 & abs(f-f_th)>=0.0005,1))]);


    tempy = [0 50 100];                     tempf = 1000*[f_th f_th f_th];           
    plot(tempf,tempy,'r');                  
    set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
    xlabel('f (mHz)');                      ylabel('Power');
end

hold off;                                   clear q tempy tempf species;

%% Plot Nullclines, State-Space trajectories
% Nullclines
R_nc = a/b;                             F_nc = c/d;     % Nullclines
 
[Rm,Fm] = meshgrid(50:10:160);     % Mesh-grid values for constructing state space
dR_num = double(subs(dR_sym,[R,F],{Rm,Fm}));
dF_num = double(subs(dF_sym,[R,F],{Rm,Fm}));
% r = ( dR_num.^2 + dF_num.^2 ).^0.5;
r=1;

figure('Name','State Space','NumberTitle','off','Position',[1 1 500 500]);          hold on;            

p_ncR = plot([0 200],[R_nc R_nc],'Color',[1 0.3 0.2],...
    'LineWidth',2,'DisplayName','N_R nullcline');
p_ncF = plot([F_nc F_nc],[0 200],'Color',[0.4 0.4 1],...
    'LineWidth',2,'DisplayName','N_F nullcline');            

% p_df = quiver(Rm,Fm,dR_num./r,dF_num./r,'g');                    % Direction field  
p_df = streakarrow(Rm,Fm,dR_num./r,dF_num./r,0.7,1);             % Direction field

p_de = plot(y_sol(:,1),y_sol(:,2),'k','LineWidth',1,...
    'DisplayName','DE trajectory');                 % plot deterministic trajectory
p_abk = plot(avg_R,avg_F,'b','LineWidth',2,...
    'DisplayName','Mean ABK trajectory');     % plot stochastic trajectory

p_ic = plot(Ri,Fi,'rp','MarkerSize',9,'MarkerFaceColor','r',...
    'DisplayName',['Initial: (' num2str(Ri) ',' num2str(Fi) ')']);     % plot initial condition
p_fp = plot(R_ss(2),F_ss(2),'oc','MarkerSize',8,...
    'MarkerFaceColor','y','DisplayName',['Center: (' num2str(R_ss(2)) ',' num2str(F_ss(2)) ')']);  % FP

set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
axis([50 160 50 160]);          set(gca,'DataAspectRatio',[1 1 1]);            
xlabel('N_R');                  ylabel('N_F');                  

leg = legend([p_ncR, p_ncF, p_fp, p_ic, p_de, p_abk]);
set(leg,'FontName','Times New Roman','FontSize',8,'Interpreter','TeX',...   
    'EdgeColor',[0.95 0.95 0.95]);                                                  hold off;

%% Confirm E is indeed a conserved quantity for deterministic solution
for v=1:size(y_sol,1)
    E(v) = double(subs(E_sym,[R F],[y_sol(v,1) y_sol(v,2)]));
end
plot(diff(E));

%% Finish
fprintf(['\n' ElapsedTime(toc)]);

clear R Ep dR* dF* R_nc F_nc;
clear h i j k m n r t v temp* tempR* texpF*;
clear fig* leg* p_*;
