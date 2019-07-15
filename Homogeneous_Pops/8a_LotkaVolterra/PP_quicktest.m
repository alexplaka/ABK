% Phase Plane of Predator-Prey system 
%     R --> 2R              Birth of R; rate constant: a (1st order)
% F + R -->  F              Death of R; rate constant: b (2nd order)
% F     -->                 Death of F; rate constant: c (1st order)
% F + R --> 2F + R          Birth of F; rate constant: d (2nd order)

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;                 clc; 
% Declare variables and functions
global a b c d;

maxTime = 100;                   % Maximum Simulation time (sec)
dt = 0.1;

agents = 600;                     % disp(['Agents = ' num2str(agents)]);
a = 1;                          % a rate constant (1st order) 
b = 0.01;                         % b rate constant (2nd order)
c = 1;                          % c rate constant (1st order) 
d = 0.01;                         % d rate constant (2nd order)

Ri = 80;                        Fi = 80;

% *** Set up differential equations of system using symbolic variables ***
syms R F positive;
dR_sym = a*R - b*F*R;
dF_sym = d*R*F - c*F;

R_ss = [0 c/d]                     
F_ss = [0 a/b]

R_nc = a/b;                             F_nc = c/d;     % Nullclines

Jac = jacobian([dR_sym,dF_sym],[R,F]);

% Nullclines
% R_nullcline_sym = solve(dR_sym == 0,R);         % No explicit solution found. Why?
% F_nullcline_sym = solve(dF_sym == 0,F);

% Linearize and calculate eigenvalues
for w=1:size(R_ss,2)   
    J = double(subs(Jac,[R,F],[R_ss(w),F_ss(w)]));
    TrJ = J(1,1) + J(2,2);
    DetJ = det(J);
    DiscrJ = TrJ^2 - 4*DetJ;
    lambda_plus(w) = (TrJ + sqrt(DiscrJ)) / 2;
    lambda_minus(w) = (TrJ - sqrt(DiscrJ)) / 2;        
end

disp([num2str(lambda_plus)]);
disp([num2str(lambda_minus)]);

% Solve differential equation
[t_sol, y_sol] = ode45(@predprey_dif,0:dt:maxTime,[Ri ; Fi]);
ul = max(max(y_sol));       % Upper Limit for plotting purposes
% Graph deterministic results
figure('Name','Time course','NumberTitle','off');
plot(t_sol,y_sol(:,1),'-b');                                      hold on;
plot(t_sol,y_sol(:,2),'-r');                               
axis tight;                             % axis([0 t_sol(end) 0 agents]);
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('time');         ylabel('# Agents');         
leg0 = legend('R deter','F deter','Location','NorthEast');        hold off;
set(leg0,'FontName','Times New Roman','FontSize',9,'Interpreter','TeX',...   
    'EdgeColor',[0.95 0.95 0.95]);                                                  hold off;

[Rm,Fm] = meshgrid(0:ul/10:ul+10);     % Mesh-grid values for constructing state space
dR_num = double(subs(dR_sym,[R,F],{Rm,Fm}));
dF_num = double(subs(dF_sym,[R,F],{Rm,Fm}));
% r = ( dR_num.^2 + dF_num.^2 ).^0.5;
r=1;

fig1 = figure('Name','State Space','NumberTitle','off','Position',[1 1 500 500]);          hold on;            

p_ncR = plot([0 ul+10],[R_nc R_nc],'Color',[1 0.3 0.2],...
    'LineWidth',2,'DisplayName','N_R nullcline');
p_ncF = plot([F_nc F_nc],[0 ul+10],'Color',[0.4 0.4 1],...
    'LineWidth',2,'DisplayName','N_F nullcline');            

% p_df = quiver(Rm,Fm,dR_num./r,dF_num./r,'g');                          
p_df = streakarrow(Rm,Fm,dR_num./r,dF_num./r,0.7,1);             % Direction field

p_de = plot(y_sol(:,1),y_sol(:,2),'k','LineWidth',1,...
    'DisplayName','DE trajectory');                 % plot deterministic trajectory
% p_abk = plot(Tr,Tf,'b');                          % plot stochastic trajectory

p_ic = plot(Ri,Fi,'rp','MarkerSize',9,'MarkerFaceColor','r',...
    'DisplayName',['Initial: (' num2str(Ri) ',' num2str(Fi) ')']);     % plot initial condition
p_fp = plot(R_ss(2),F_ss(2),'oc','MarkerSize',8,...
    'MarkerFaceColor','y','DisplayName',['Center: (' num2str(R_ss(2)) ',' num2str(F_ss(2)) ')']);  % FP

set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
axis([0 ul+10 0 ul+10]);          set(gca,'DataAspectRatio',[1 1 1]);            
xlabel('N_R');                  ylabel('N_F');                  

leg1 = legend([p_ncR, p_ncF, p_fp, p_ic, p_de]);
set(leg1,'FontName','Times New Roman','FontSize',9,'Interpreter','TeX',...   
    'EdgeColor',[0.95 0.95 0.95],'Location','Best');
hold off;

%% Spectral Analysis
fs = 1/dt;                                          % sampling frequency
t_osc = maxTime/10:dt:maxTime;                      
R_osc = y_sol(t_osc(1)/dt+1:end,1)';
F_osc = y_sol(t_osc(1)/dt+1:end,2)';

Osc = [R_osc; F_osc];
f_th = abs(imag(lambda_plus(2))) / (2*pi);             % Theoretical oscillation frequency (Hz)
period_th = 1 / f_th;                               % Theoretical period (sec)

figure('Name','Spectral Analysis: DE','NumberTitle','off','Position',[1 1 350 800]);          

% do FFT on R, F data
n = 2^(nextpow2(length(R_osc))+4);                      
f = fs / n * (0:n/2);       

species = {'\underbar{Species R}','\underbar{Species F}'}; 

Z = zeros(size(species,2),n);               Zn = zeros(size(species,2),n);

for b=1:size(species,2)
    Z(b,:) = fft(Osc(b,:),n);
    Zn(b,:) = abs(Z(b,:)/n);
    subplot(size(species,2),1,b);                                             hold on;
    title(species{b},'FontName','Times New Roman','FontSize',12,'Interpreter','LaTex');  
    plot(1000*f,Zn(b,1:n/2+1));                                   
    
%     axis([0.5*f_th 2*f_th 0 12*b]);
    xlim([500*f_th 1500*f_th]);
    ylim([0 1.5*Zn(b,find(abs(f-f_th)<=0.0002 & abs(f-f_th)>=0.00005,1))]);


    tempy = [0 50 100];                     tempf = 1000*[f_th f_th f_th];           
    plot(tempf,tempy,'r');                  
    set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
    xlabel('f (mHz)');                      ylabel('Power');
end
  
hold off;                               clear b tempy tempf;