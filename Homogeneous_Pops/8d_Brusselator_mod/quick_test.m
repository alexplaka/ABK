% For testing ODE of the modified Brusselator system.

%   --> X           (birth of X - 0th order; rate constant: b)
% X -->             (death or degradadation of X - 1st order wrt X; rate constant: d)
% X --> Y           (1st order wrt X; rate constant: c)

% Replace the 3rd order reaction
% 2X + Y --> 3X     (3rd order: 2nd wrt X, 1st wrt Y; rate constant: a)
% with: 
% 2X <--> Z                     k_f/k_r = K_eq
% Z + Y --> Z + X               2nd order rate constant a 

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;  clc;

global a b c d k_f k_r;

maxTime = 5000;                         % Maximum Simulation time (sec)
dt = 1;                                 

a = 0.01;                               % Z + Y --> Z + X rate constant (2nd order)
b = 0.3;                               % birth rate constant (0th order) for X
c = 0.075;                                % X --> Y conversion rate constant (1st order)
d = 0.02;                               % death rate constant (1st order) for X
k_f = 0.01;
k_r = 1;
K = k_f / k_r;

Xi = 0;         Yi = 0;         Zi = 0;         % ** Initial condition **

syms X Y Z positive;

DEchoice = input('Canonical (0) or Agent-Based DEs (1):');
% DEchoice = 1;                             % For testing

if DEchoice == 0   
    [t_sol, y_sol] = ode45(@bruss_can_dif,0:maxTime/1000:maxTime,[Xi ; Yi ; Zi]);   % Canonical DEs
    
    dX_sym = b - 2 * k_f * X^2 + 2 * k_r * Z + a * Z * Y  - (c + d) * X; 
    dY_sym = c * X - a * Z * Y;
    dZ_sym = k_f * X^2 - k_r * Z;
    
    X_nc_sym = ((c + d)*X - b) / (a*K * X^2);   % X nullcline. Given in form X = expression
    Y_nc_sym = c / (a*K*X);                     % Y nullcline. Given in form Y = expression
    Z_nc_sym = K * X^2;                         % Z nullcline. Given in form Z = expression
    
    X_ss = b / d;                  Y_ss = (c * d) / (b * a*K);         Z_ss = K * b^2 / d^2;
        
    startvalue = 1;                             % For making mesh-grid
    
elseif DEchoice == 1
    [t_sol, y_sol] = ode45(@bruss_ab_dif,0:maxTime/1000:maxTime,[Xi ; Yi ; Zi]);    % Agent-Based DEs
    
    dX_sym = b - 2 * k_f * X * (X-1) + 2 * k_r * Z + a * Z * Y - (c + d) * X;
    dY_sym =  c * X - a * Z * Y;
    dZ_sym = k_f * X * (X-1) - k_r * Z;
    
    X_nc_sym = ((c + d)*X - b) / (a*K * X*(X-1));   % X nullcline. Given in form X = expression
    Y_nc_sym = c / (a*K * (X-1));                   % Y nullcline. Given in form Y = expression
    Z_nc_sym = K * X * (X-1);                       % Z nullcline. Given in form Z = expression
    
    X_ss = b / d;              Y_ss = (c * d) / (a*K * (b - d));         Z_ss = K * b^2 / d^2;
        
    startvalue = 2;                             % For making mesh-grid; X ~= 1 in Agent-Based DE

else
    disp('Wrong choice!');
end

disp(['Theoretical steady state X = ' num2str(X_ss) ...
    ' , Y = ' num2str(Y_ss) ' , Z = ' num2str(Z_ss)]);

% Linearize and calculate eigenvalues
Jac = jacobian([dX_sym,dY_sym,dZ_sym],[X,Y,Z]);
J = eval(subs(Jac,[X,Y,Z],[X_ss,Y_ss,Z_ss]));          
lambda = eig(J);

disp(['lambda = ' num2str(lambda')]);

if imag(lambda(2)) ~= 0 
    f_th = abs(imag(lambda(2))) / (2*pi);             % Theoretical oscillation frequency (Hz)
    period_th = 1 / f_th;                             % Theoretical period (sec)
    disp(['Frequency = ' num2str(f_th) ' Hz']);
    disp(['Period    = ' num2str(period_th) ' sec']);
end

%% Plot DE time course
fig0 = figure('Name','DE Time Course ',...
    'NumberTitle','off','Position',[1 1 500 500]);                          hold on;
title(['N_{X,i} = ' num2str(Xi) ' , N_{Y,i} = ' num2str(Yi) ' , N_{Z,i} = ' num2str(Zi)],...
    'FontSize',12,'FontName','Times New Roman');

p_Xd = plot(t_sol,y_sol(:,1),'--','Color',[0.75 0 1],'LineWidth',2,'DisplayName','DE N_{X}(t)');                                       
p_Yd = plot(t_sol,y_sol(:,2),'--','Color',[0 0.5 0],'LineWidth',2,'DisplayName','DE N_{Y}(t)');
p_Zd = plot(t_sol,y_sol(:,3),'--','Color',[0 0 0.75],'LineWidth',2,'DisplayName','DE N_{Z}(t)');

axis tight;                                     
% axis([0 time(end) 0 50]);                               
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');                           hold off;
xlabel('t (sec)');                 ylabel('N(t)'); 

leg0 = legend([p_Xd , p_Yd , p_Zd]);           
set(leg0,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');       

%% Calculate Nullclines, meshgrid for direction field
cap = max(max(y_sol))*1.2;                          % ** Set appropriate state-space window **

% X_nc = double(subs(X_nc_sym,X,startvalue:cap));      
% Y_nc = double(subs(Y_nc_sym,X,startvalue:cap));
% Z_nc = double(subs(Z_nc_sym,X,startvalue:cap));

[Xm,Ym,Zm] = meshgrid(startvalue:cap/10:cap);       % Mesh-grid values for constructing state space

dX_num = double(subs(dX_sym,[X,Y,Z],{Xm,Ym,Zm}));
dY_num = double(subs(dY_sym,[X,Y,Z],{Xm,Ym,Zm}));
dZ_num = double(subs(dZ_sym,[X,Y,Z],{Xm,Ym,Zm}));

% r = ( dX_num.^2 + dY_num.^2 + dZ_num.^2 ).^0.5;
r=1;
%% Plot Nullclines, State-Space trajectories
figure('Name','State Space','NumberTitle','off','Position',[1 1 500 500]);          hold on;            

% p_ncX = plot(startvalue:cap,X_nc,'Color',[1 0.3 0.2],...
%     'LineWidth',2,'DisplayName','N_X nullcline');
% p_ncY = plot(startvalue:cap,Y_nc,'Color',[0.4 0.4 1],...
%     'LineWidth',2,'DisplayName','N_Y nullcline');            
% p_ncZ = plot(startvalue:cap,Z_nc,'Color',[0.4 1 0.3],...
%     'LineWidth',2,'DisplayName','N_Z nullcline');            

quiver3(Xm,Ym,Zm,dX_num./r,dY_num./r,dZ_num./r,'g');      % plot direction field                           
% streakarrow(Xm,Ym,dX_num./r,dY_num./r,0.7,1);           % plot direction field (2D only)

p_de = plot3(y_sol(:,1),y_sol(:,2),y_sol(:,3),'k','LineWidth',1,...
    'DisplayName','DE trajectory');                       % plot deterministic trajectory

p_ic = plot3(Xi,Yi,Zi,'rp','MarkerSize',9,'MarkerFaceColor','r',...  % plot initial condition
    'DisplayName',['Initial: (' num2str(Xi) ',' num2str(Yi) ',' num2str(Zi) ')']);
p_fp = plot(X_ss,Y_ss,'oc','MarkerSize',8,...
    'MarkerFaceColor','y','DisplayName',...
    ['FP: (' num2str(X_ss) ',' num2str(Y_ss) ',' num2str(Z_ss) ')']);             % FP

set(gca,'XMinorTick','on','YMinorTick','on','ZMinorTick','on','Box','off');
axis([0 cap 0 cap 0 cap]);            set(gca,'DataAspectRatio',[1 1 1]);            
xlabel('N_X');                  ylabel('N_Y');                  zlabel('N_Z');

% leg = legend([p_ncX, p_ncY, p_fp, p_ic, p_de]);
leg = legend([p_fp, p_ic, p_de]);
set(leg,'FontName','Times New Roman','FontSize',8,'Interpreter','TeX',...   
    'EdgeColor',[0.95 0.95 0.95]);                                                  hold off;
