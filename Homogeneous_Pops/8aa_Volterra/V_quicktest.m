% Predator-Prey system (Volterra system)
%     R --> 2R              Density-dependent Birth of R; rate constant: a (1st order)
% F + R -->  F              Death of R; rate constant: b (2nd order)
% F     -->                 Death of F; rate constant: c (1st order)
% F + R --> 2F + R          Birth of F; rate constant: d (2nd order)

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;                 clc; 
% Declare variables and functions
global a b c d K;

maxTime = 100;                   % Maximum Simulation time (sec)
dt = 0.1;

agents = 200;                     % disp(['Agents = ' num2str(agents)]);

a = 1;                          % a rate constant (1st order) 
b = 0.01;                         % b rate constant (2nd order)
c = 1;                          % c rate constant (1st order) 
d = 0.01;                         % d rate constant (2nd order)
K = 120;                        % Carrying capacity for R

Ri = 100;                        Fi = 100;

% *** Set up differential equations of system using symbolic variables ***
syms R F positive;
dR_sym = a*R*(1-R/K) - b*F*R;
dF_sym = d*R*F - c*F;

R_ss = [0 c/d]                     
F_ss = [0 a/b * (1 - c/(d*K))]

%% Stability analysis
% ** Nullclines ** (explicitly stated)
R_nc = a/b .* (1 - [0:agents]/K);   % a/K * R^2 + R*(bF-a) = 0  implies F = a/b * (1 - R/K)               
F_nc = c/d;                         % R = c/d;

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

disp([num2str(lambda_plus)]);
disp([num2str(lambda_minus)]);

% Conditions for having a Stable Spiral
cond1 = 4*d*K * (d*K/c - 1);
if a < cond1,    disp('Stable Spiral predicted (wrt param a).');   end

cond2 = ( 2*d/a * (sqrt(1 + a/c) - 1) )^-1;
if K > cond2,    disp('Stable Spiral predicted (wrt param K).');   end

%% Solve differential equation and plot
[t_sol, y_sol] = ode45(@predprey_dif,0:dt:maxTime,[Ri ; Fi]);

% ul = max(max(y_sol));       % Upper Limit for plotting purposes
ul = agents;

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
%% Phase Plane
% [Rm,Fm] = meshgrid(ul/4:ul/15:3*ul/4);          % Mesh-grid values for constructing state space   
% 
% dR_num = double(subs(dR_sym,[R,F],{Rm,Fm}));
% dF_num = double(subs(dF_sym,[R,F],{Rm,Fm}));
% % r = ( dR_num.^2 + dF_num.^2 ).^0.5;
% r=1;
% 
% fig1 = figure('Name','State Space','NumberTitle','off','Position',[1 1 500 500]);          hold on;            
% 
% p_ncR = plot(0:agents,R_nc,'Color',[1 0.3 0.2],...
%     'LineWidth',2,'DisplayName','N_R nullcline');
% p_ncF = plot([F_nc F_nc],[0 ul],'Color',[0.4 0.4 1],...
%     'LineWidth',2,'DisplayName','N_F nullcline');            
% 
% % p_df = quiver(Rm,Fm,dR_num./r,dF_num./r,'g');                          
% p_df = streakarrow(Rm,Fm,dR_num./r,dF_num./r,0.7,1);             % Direction field
% 
% p_de = plot(y_sol(:,1),y_sol(:,2),'k','LineWidth',1,...
%     'DisplayName','DE trajectory');                 % plot deterministic trajectory
% % p_abk = plot(Tr,Tf,'b');                          % plot stochastic trajectory
% 
% p_ic = plot(Ri,Fi,'rp','MarkerSize',9,'MarkerFaceColor','r',...
%     'DisplayName',['Initial: (' num2str(Ri) ',' num2str(Fi) ')']);     % plot initial condition
% p_fp = plot(R_ss(2),F_ss(2),'oc','MarkerSize',8,...
%     'MarkerFaceColor','y',...
%     'DisplayName',['Stable Spiral: (' num2str(R_ss(2)) ',' num2str(F_ss(2)) ')']);  % FP
% 
% set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
% axis([ul/4 3*ul/4 ul/4 3*ul/4]);              set(gca,'DataAspectRatio',[1 1 1]);            
% xlabel('N_R');                  ylabel('N_F');                  
% 
% leg1 = legend([p_ncR, p_ncF, p_fp, p_ic, p_de]);
% set(leg1,'FontName','Times New Roman','FontSize',9,'Interpreter','TeX',...   
%     'EdgeColor',[0.95 0.95 0.95]);                                                  hold off;
