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

clear;                   tic;
clc; 
rng(0);

% Declare variables and functions
global a b c d K;

maxTime = 100;                     % Maximum Simulation time (sec)
dt = 0.01;                         % Constant (fixed) time step increment (sec)
t_steps = maxTime / dt;

agents = 500;                     % Max number of agents allowed
% If a population size exceeds this, then species F is extinct and R grows exponentially.

a = 1;                          % a rate constant 
b = 0.01;                       % b rate constant (2nd order)
c = 1;                          % c rate constant (1st order) 
d = 0.01;                       % d rate constant (2nd order)
K = agents;                     % Carrying capacity for R

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

%% Initialize - Preallocate memory for variables; ** Initial conditions **
Pb = zeros(1,t_steps);              % For storing probability value at each time step
Pd = zeros(1,t_steps);              % For storing probability value at each time step

Tr = zeros(1,t_steps);              % Total sum of R agents in each time step	
Tf = zeros(1,t_steps);              % Total sum of F agents in each time step
time = zeros(1,t_steps);

% ******** Initial conditions - Number of R, Rp Agents ********
Tr(1) = Ri;                		Tf(1) = Fi;					
% ************************************************************

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

%% ABK Simulation
t = 1;                                      % Time counter variable

Pc = 1 - exp(-c * dt);                      % Pr for death of F: F --> 

while t*dt <= maxTime && Tr(t) < 0.95*K
    
    Pa(t) = a * (1 - Tr(t)/K) * dt;         % Pr for birth of R: R --> 2R
    Pb(t) = b * Tf(t) * dt;                 % Pr for death of R: F + R --> F  (wrt R)
    Pd(t) = d * Tr(t) * dt;                 % Pr for birth of F: F + R --> 2F + R  (wrt F)
        
    tempR_al = find(Rv(1,:)==1);            tempR_dd = find(Rv(1,:)==0);
    
    for i = 1:size(tempR_al,2)
        r = rand;
        if r < Pb(t)
            Rv(2,tempR_al(i)) = 0;         % R(i) agent dies     (2nd order)         
        elseif r >= Pb(t) && r < Pb(t)+Pa(t)   % R agent 'i' gives birth 
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
    
    Tr(t+1) = sum(Rv(2,:));                  Tf(t+1) = sum(Fv(2,:));
    time(t+1) = time(t) + dt;
    
    Rv(1,:) = Rv(2,:);                       Fv(1,:) = Fv(2,:); 
    t = t + 1;  
    
end

% Remove unnecessary terminal 0's from arrays
if t < t_steps
    Tr = Tr(1:t-1);                         Tf = Tf(1:t-1);   
    time = time(1:t-1);
end

% disp(['ABK sim terminal R value  = ' num2str(Tr(end))]);  
% disp(['ABK sim terminal F value  = ' num2str(Tf(end))]);        
clear Rv Fv;
finaltime = (t-1) * dt;

% Upper Limit for plotting purposes
% ul = max(Tr);                       % ul = max(max(y_sol));       
ul = agents;

%% Solve differential equation (was waiting for variable 'finaltime')
[t_sol, y_sol] = ode45(@predprey_dif,0:finaltime/1000:finaltime,[Ri ; Fi]);
% disp(['Diff eq terminal R value  = ' num2str(y_sol(end,1))]);
% disp(['Diff eq terminal F value  = ' num2str(y_sol(end,2))]);

%% Graph ABK (stochastic) and deterministic results
fig0 = figure('Name','Time course','NumberTitle','off','Position',[1 1 500 500]);      hold on;

p_Rs = plot(time,Tr,'r','DisplayName','ABK N_{R}(t)');                                        
p_Fs = plot(time,Tf,'g','DisplayName','ABK N_{F}(t)');
p_Rd = plot(t_sol,y_sol(:,1),'--','Color',[0.75 0 1],'LineWidth',2,...
    'DisplayName','DE N_{R}(t)');                                       
p_Fd = plot(t_sol,y_sol(:,2),'--','Color',[0 0.5 0],'LineWidth',2,...
    'DisplayName','DE N_{F}(t)');

% axis tight;                         
axis([0 t_sol(end) 0 200]);          
xlabel('t (months)');         ylabel('N(t)');
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');                       hold off;

leg0 = legend([p_Rs , p_Fs , p_Rd , p_Fd]);           % All species
% leg0 = legend([p_Rd , p_Fd]);
set(leg0,'FontName','Times New Roman','FontSize',8,'Interpreter','TeX',...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');                                              

% Create textbox
annotation(fig0,'textbox',[0.03 0.92 0.05 0.07],'String',{'b)'},'FontSize',12,...
    'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');

%% Plot State Space: nullclines, direction field, trajectories
% ** Nullclines ** (explicitly stated)
R_nc = a/b .* (1 - [0:5:agents]/K);   % a/K * R^2 + R*(bF-a) = 0  implies F = a/b * (1 - R/K)               
F_nc = c/d;                         % R = c/d;

[Rm,Fm] = meshgrid(ul/4:ul/15:3*ul/4);        % Mesh-grid values for constructing state space
dR_num = double(subs(dR_sym,[R,F],{Rm,Fm}));
dF_num = double(subs(dF_sym,[R,F],{Rm,Fm}));
% r = ( dR_num.^2 + dF_num.^2 ).^0.5;
r=1;
%
fig1 = figure('Name','State Space','NumberTitle','off','Position',[1 1 500 500]);          hold on;            

p_ncR = plot(0:5:agents,R_nc,'Color',[1 0.3 0.2],...
    'LineWidth',2,'DisplayName','N_R nullcline');
p_ncF = plot([F_nc F_nc],[0 ul],'Color',[0.4 0.4 1],...
    'LineWidth',2,'DisplayName','N_F nullcline');            

% p_df = quiver(Rm,Fm,dR_num./r,dF_num./r,'g');                          
p_df = streakarrow(Rm,Fm,dR_num./r,dF_num./r,0.7,1);             % Direction field

p_de = plot(y_sol(:,1),y_sol(:,2),'k','LineWidth',1,...
    'DisplayName','DE trajectory');         % plot deterministic trajectory
p_abk = plot(Tr,Tf,'b','LineWidth',1,...
    'DisplayName','DE trajectory');         % plot stochastic trajectory

p_ic = plot(Ri,Fi,'rp','MarkerSize',9,'MarkerFaceColor','r',...
    'DisplayName',['Initial: (' num2str(Ri) ',' num2str(Fi) ')']);     % plot initial condition
p_fp = plot(R_ss(2),F_ss(2),'oc','MarkerSize',8,...
    'MarkerFaceColor','y','DisplayName',['Center: (' num2str(R_ss(2)) ',' num2str(F_ss(2)) ')']);  % FP

set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
% axis([0 ul 0 ul]);          
axis([ul/10 3*ul/10 ul/10 3*ul/10]);
set(gca,'DataAspectRatio',[1 1 1]);            
xlabel('N_R');                ylabel('N_F');                  

leg1 = legend([p_ncR, p_ncF, p_fp, p_ic, p_abk, p_de]);
set(leg1,'FontName','Times New Roman','FontSize',9,'Interpreter','TeX',...   
    'EdgeColor',[0.95 0.95 0.95]);                                                  hold off;

% Create textbox
annotation(fig1,'textbox',[0.03 0.92 0.05 0.07],'String',{'a)'},'FontSize',12,...
    'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');

%% Spectral Analysis for ABK trajectory
fs = 1/dt;                                          % sampling frequency
t_osc = maxTime/10:dt:maxTime;                      
R_osc = Tr(1,t_osc(1)/dt+1:end);            F_osc = Tf(1,t_osc(1)/dt+1:end);

Osc = [R_osc; F_osc];
f_th = abs(imag(lambda_plus(2))) / (2*pi);          % Theoretical oscillation frequency (Hz)
period_th = 1 / f_th;                               % Theoretical period (sec)

figure('Name','Spectral Analysis: Trial','NumberTitle','off','Position',[1 1 350 800]);          

% do FFT on R, F data
n = 2^(nextpow2(length(R_osc))+4);                      
f = fs / n * (0:n/2);       

species = {'\underbar{Species R}','\underbar{Species F}'}; 

Z = zeros(size(species,2),n);               Zn = zeros(size(species,2),n);

for q=1:size(species,2)
    Z(q,:) = fft(Osc(q,:),n);
    Zn(q,:) = abs(Z(q,:)/n);
    subplot(size(species,2),1,q);                                             hold on;
    title(species{q},'FontName','Times New Roman','FontSize',12,'Interpreter','LaTex');  
    plot(1000*f,Zn(q,1:n/2+1));                                   
    
%     axis([0.5*f_th 2*f_th 0 12*q]);
    xlim([500*f_th 1500*f_th]);
    ylim([0 5*Zn(q,find(abs(f-f_th)<=0.002 & abs(f-f_th)>=0.0005,1))]);


    tempy = [0 50 100];                     tempf = 1000*[f_th f_th f_th];           
    plot(tempf,tempy,'r');                  
    set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
    xlabel('f (mHz)');                      ylabel('Power');
end
  
hold off;                               clear q tempy tempf;
%% Finish
clear r w i j k;
toc
