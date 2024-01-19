% For testing ODE of the Brusselator system.

%   --> X           (birth of X - 0th order; rate constant: b)
% X -->             (death or degradadation of X - 1st order wrt X; rate constant: d)
% X --> Y           (1st order wrt X; rate constant: c)
% 2X + Y --> 3X     (3rd order: 2nd wrt X, 1st wrt Y; rate constant: a)

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear; % clc;

global a b c d;

maxTime = 2000;                          % Maximum Simulation time (sec)
dt = 1;                                 

a = 0.0005;                           % 2 X + Y --> 3 X rate constant (3rd order)
b = 0.5;                              % birth rate constant (0th order) for X
c = 0.11;                             % X --> Y conversion rate constant (1st order)
d = 0.05;                             % death rate constant (1st order) for X

Xi = 0;         Yi = 0;                % ** Initial condition **

syms X Y positive;

DEchoice = input('Canonical (0) or Agent-Based DEs (1):');

if DEchoice == 0   
    [t_sol, y_sol] = ode45(@bruss_can_dif,0:maxTime/1000:maxTime,[Xi ; Yi]);   % Canonical DEs
    
    dX_sym = b + a * X^2 * Y - (c + d) * X;
    dY_sym =  c * X - a * X^2 * Y;
    
    X_nc_sym = ((c + d)*X - b) / (a*X^2);       % X nullcline. Given in form Y = expression
    Y_nc_sym = c / (a*X);                       % Y nullcline. Given in form Y = expression

    X_ss = b / d;                        Y_ss = (c * d) / (b * a); 
    
    % Condition for existence of limit cycle (Canonical DEs):
    crit = d + a * b^2 / d^2;
    if c > crit,    disp('Limit cycle expected (Caconical DEs)');      end
    
    startvalue = 1;                             % For making mesh-grid
    
elseif DEchoice == 1
    [t_sol, y_sol] = ode45(@bruss_ab_dif,0:maxTime/1000:maxTime,[Xi ; Yi]);    % Agent-Based DEs
    
    dX_sym = b + a * X * (X-1) * Y - (c + d) * X;
    dY_sym =  c * X - a * X * (X-1) * Y;
    
    X_nc_sym = ((c + d)*X - b) / (a*X*(X-1));   % X nullcline. Given in form Y = expression
    Y_nc_sym = c / (a*(X-1));                   % Y nullcline. Given in form Y = expression

    X_ss = b / d;                       Y_ss = (c * d) / (a*(b - d));
    
    % Condition for existence of limit cycle (Agent-based DEs)
    crit = ( d + a * b * (b-d) / d^2 ) / ( b / (b-d) );
    if c > crit,    disp('Limit cycle expected (Agent-Based DEs)');    end
    
    startvalue = 2;                             % For making mesh-grid; X ~= 1 in Agent-Based DE

else
    disp('Wrong choice!');
end

disp(['Theoretical steady state X = ' num2str(X_ss) ' , Y = ' num2str(Y_ss)]);

% Linearize and calculate eigenvalues
Jac = jacobian([dX_sym,dY_sym],[X,Y]);
J = eval(subs(Jac,[X,Y],[X_ss,Y_ss]));          
TrJ = J(1,1) + J(2,2);
DetJ = det(J);
DiscrJ = TrJ^2 - 4*DetJ;

lambda(1) = (TrJ + sqrt(DiscrJ)) / 2;
lambda(2) = (TrJ - sqrt(DiscrJ)) / 2;        
disp(['lambda = ' num2str(lambda(1))]);

if real(lambda(1))>0 && imag(lambda(1)) ~= 0 
    f_th = abs(imag(lambda(1))) / (2*pi);             % Theoretical oscillation frequency (Hz)
    period_th = 1 / f_th;                             % Theoretical period (sec)
    disp(['Frequency = ' num2str(f_th) ' Hz']);
    disp(['Period    = ' num2str(period_th) ' sec']);
end

%% Plot DE time course
fig0 = figure('Name','DE Time Course ',...
    'NumberTitle','off','Position',[1 1 500 500]);                          hold on;
title(['N_{X,i} = ' num2str(Xi) ' , N_{Y,i} = ' num2str(Yi)],...
    'FontSize',12,'FontName','Times New Roman');

p_Xd = plot(t_sol,y_sol(:,1),'--','Color',[0.75 0 1],'LineWidth',2,'DisplayName','DE N_{X}(t)');                                       
p_Yd = plot(t_sol,y_sol(:,2),'--','Color',[0 0.5 0],'LineWidth',2,'DisplayName','DE N_{Y}(t)');

axis tight;                                     
% axis([0 time(end) 0 50]);                               
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');                           hold off;
xlabel('t (sec)');                 ylabel('N(t)'); 

leg0 = legend([p_Xd , p_Yd]);           
set(leg0,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');       

%% Calculate Nullclines, meshgrid for direction field
cap = max(max(y_sol))*1.2;                 % ** Set appropriate state-space window **

X_nc = double(subs(X_nc_sym,X,startvalue:cap));      
Y_nc = double(subs(Y_nc_sym,X,startvalue:cap));

[Xm,Ym] = meshgrid(startvalue:cap/10:cap);     % Mesh-grid values for constructing state space

dX_num = double(subs(dX_sym,[X,Y],{Xm,Ym}));
dY_num = double(subs(dY_sym,[X,Y],{Xm,Ym}));
% r = ( dX_num.^2 + dY_num.^2 ).^0.5;
r=1;
%% Plot Nullclines, State-Space trajectories
figure('Name','State Space','NumberTitle','off','Position',[1 1 500 500]);          hold on;            

p_ncX = plot(startvalue:cap,X_nc,'Color',[1 0.3 0.2],...
    'LineWidth',2,'DisplayName','N_X nullcline');
p_ncY = plot(startvalue:cap,Y_nc,'Color',[0.4 0.4 1],...
    'LineWidth',2,'DisplayName','N_Y nullcline');            

% quiver(Xm,Ym,dX_num./r,dY_num./r,'g');                  % plot direction field                           
streakarrow(Xm,Ym,dX_num./r,dY_num./r,0.7,1);           % plot direction field

p_de = plot(y_sol(:,1),y_sol(:,2),'k','LineWidth',1,...
    'DisplayName','DE trajectory');                 % plot deterministic trajectory

p_ic = plot(Xi,Yi,'rp','MarkerSize',9,'MarkerFaceColor','r',...
    'DisplayName',['Initial: (' num2str(Xi) ',' num2str(Yi) ')']);     % plot initial condition
p_fp = plot(X_ss,Y_ss,'oc','MarkerSize',8,...
    'MarkerFaceColor','y','DisplayName',['FP: (' num2str(X_ss) ',' num2str(Y_ss) ')']);  % FP

set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
axis([0 cap 0 cap]);            set(gca,'DataAspectRatio',[1 1 1]);            
xlabel('N_X');                  ylabel('N_Y');                  

leg = legend([p_ncX, p_ncY, p_fp, p_ic, p_de]);
set(leg,'FontName','Times New Roman','FontSize',8,'Interpreter','TeX',...   
    'EdgeColor',[0.95 0.95 0.95]);                                                  hold off;
