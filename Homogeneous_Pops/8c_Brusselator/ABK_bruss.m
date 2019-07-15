% The Brusselator (a hypothetical system) reaction scheme (reference?). 
% ABK implementation.

%        ^                               b  = birth (0th order) rate constant for X
%      d |                               d  = death (1st order) rate constant for X
%        |     c                         c  = X --> Y (1st order) rate constant 
%   -->  X   <-->  Y                     a =  2 X + Y --> 3 X rate constant (3rd order)
%    b         a                        
% 
% Or, alternatively stated: 
% 
%   --> X           (birth of X - 0th order; rate constant: b)
% X -->             (death or degradadation of X - 1st order wrt X; rate constant: d)
% X --> Y           (1st order wrt X; rate constant: c)
% 2X + Y --> 3X     (3rd order: 2nd wrt X, 1st wrt Y; rate constant: a)

% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;          clc;         tic;
 
% Declare variables and functions
global a b c d;

maxTime = 50000;                          % Maximum Simulation time (sec)
dt = 0.02;                              % Constant (fixed) time step increment (sec)
t_steps = maxTime / dt;

a = 0.00003;                           % 2 X + Y --> 3 X rate constant (3rd order)
b = 0.2;                              % birth rate constant (0th order) for X
c = 0.02;                             % X --> Y conversion rate constant (1st order)
d = 0.01;                             % death rate constant (1st order) for X

Xi = 0;         Yi = 0;                % ** Initial condition **

agents = 200;                       % disp(['Agents = ' num2str(agents)]);

%% Initialize - Preallocate memory for variables; ** Initial conditions **
Sx = zeros(1,t_steps);              % Sum of X agents in each time step	
Sy = zeros(1,t_steps);              % Sum of Y agents in each time step

% ******** Initial conditions - Number of A, B agents ********
Sx(1) = Xi;                         Sy(1) = Yi;					
% ************************************************************

tempX = zeros(1,agents);                        
% Put "1" where agents are "alive", then randomize the array
for w=1:Sx(1),                      tempX(w)=1;             end
tempX = RandArray(tempX);           % Randomize array
Xv = [tempX ; tempX];               % Initialize vector storing state of X agents     

tempY = zeros(1,agents); 
% Put "1" where agents are "alive", then randomize the array
for z=1:Sy(1),                      tempY(z)=1;             end
tempY = RandArray(tempY);           % Randomize array
Yv = [tempY ; tempY];               % Initialize vector storing state of B agents
% Notes on Xv, Yv:
% - Markov process, so only previous and current time steps needed --> 2 rows:          
clear w z tempX tempY;

%% ABK Simulation
% rxns = zeros(1,t_steps);
P_a = zeros(1,t_steps);          % Preallocate memory

% Time-independent probabilities
P_b = b * dt;                    % Probability of "birth of X" process (0th order)
P_d = 1 - exp(-d * dt);          % Probability of "death" process (1st order for X)
P_c = 1 - exp(-c * dt);          % Probability of 1st order rx X --> Y (1st order wrt X)

for t=1:t_steps
        
    % Take care of 0th order process first
    if rand < P_b                                   % "Birth", 0th order reaction
        q = find(Xv(1,:)==0);                       % Find X agents who are in "dead" state 
        Xv(2,q(ceil(rand * size(q,2)))) = 1;        % an X agent is born
    end                                    
    
    % Degradation of X    or    X --> Y
    temp1 = find(Xv(1,:)==1);
    for i=1:size(temp1,2)
        r = rand;
        if r < P_d
            Xv(2,temp1(i)) = 0;                     % X agent dies
        elseif r >= P_d && r < P_d+P_c
            Xv(2,temp1(i)) = 0;                     % X agent dies
            temp2 = find(Yv(1,:)==0);               % Find Y agents who are in "dead" state
            p = ceil(rand*size(temp2,2));           % Pick Y "dead" agent randomly
            Yv(2,temp2(p)) = 1;                     % Y agent is born
        end
    end
    
    % 2 X + Y --> 3 X , evaluate P wrt Y
    P_a(t) = a * Sx(t) * (Sx(t)-1) * dt;          % Prob of 2 X + Y --> 3 X          
    temp3 = find(Yv(1,:)==1);
    for j=1:size(temp3,2)
%         Xalive = Sx(t) - 2*rxns(t);  
%         P_a(t) = a * Xalive * (Xalive-1) * dt;      % Prob of 2 X + Y --> 3 X          
        if rand < P_a(t)
            Yv(2,temp3(j)) = 0;                     % Y agent dies
            temp4 = find(Xv(1,:)==0);               % Find X agents who are "dead"
            q = ceil(rand*size(temp4,2));           % Pick X "dead" agent randomly
            Xv(2,temp4(q)) = 1;                     % X agent is born 
%             rxns(t) = rxns(t) + 1;
        end
    end
            
    Sx(t+1) = sum(Xv(2,:));                  Sy(t+1) = sum(Yv(2,:));
    
    Xv(1,:) = Xv(2,:);                       Yv(1,:) = Yv(2,:); 
    
end

%% Symbolic calculations and numerical solution to DEs
syms X Y positive;

DEchoice = input('Canonical (0) or Agent-Based DEs (1):');

if DEchoice == 0   
    [t_sol, y_sol] = ode45(@bruss_can_dif,0:maxTime/1000:maxTime,[Xi ; Yi]);   % Canonical DEs
    
    dX_sym = b + a * X^2 * Y - (c + d) * X;
    dY_sym =  c * X - a * X^2 * Y;
    
    X_nc_sym = ((c + d)*X - b) / (a*X^2);       % X nullcline. Given in form Y = expression
    Y_nc_sym = c / (a*X);                       % Y nullcline. Given in form Y = expression

    X_ss = b / d;                          Y_ss = (c * d) / (b * a);  
    
    % Condition for existence of limit cycle (Canonical DEs):
    crit = d + a * b^2 / d^2;
    if c > crit,    disp('Limit cycle expected (caconical DEs)');      end
    
    startvalue = 1;                             % For making mesh-grid
    
elseif DEchoice == 1
    [t_sol, y_sol] = ode45(@bruss_ab_dif,0:maxTime/1000:maxTime,[Xi ; Yi]);    % Agent-Based DEs
    
    dX_sym = b + a * X * (X-1) * Y - (c + d) * X;
    dY_sym =  c * X - a * X * (X-1) * Y;
    
    X_nc_sym = ((c + d)*X - b) / (a*X*(X-1));   % X nullcline. Given in form Y = expression
    Y_nc_sym = c / (a*(X-1));                   % Y nullcline. Given in form Y = expression
    
    X_ss = b / d;                          Y_ss = (c * d) / (a*(b - d));
    
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
TrJ = J(1,1) + J(2,2);                  DetJ = det(J);
DiscrJ = TrJ^2 - 4*DetJ;

lambda(1) = (TrJ + sqrt(DiscrJ)) / 2;
lambda(2) = (TrJ - sqrt(DiscrJ)) / 2;        
disp(['lambda = ' num2str(lambda)]);

if real(lambda(1))>0 && imag(lambda(1)) ~= 0 
    f_th = abs(imag(lambda(1))) / (2*pi);             % Theoretical oscillation frequency (Hz)
    period_th = 1 / f_th;                             % Theoretical period (sec)
    disp(['Frequency = ' num2str(f_th) ' Hz']);
    disp(['Period    = ' num2str(period_th) ' sec']);
end

%% Graph ABK (stochastic) and deterministic results
figure('Name','Sample Time course','NumberTitle','off');
time = 0:dt:maxTime;                                hold on;
plot(time,Sx,'b');                                        
plot(time,Sy,'r');
plot(t_sol,y_sol(:,1),'--c','LineWidth',3);                  
plot(t_sol,y_sol(:,2),'--m','LineWidth',3);                               

% axis([0 t_sol(end) 0 agents]);          
axis tight;
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');           hold off;

xlabel('time');         ylabel('# Agents');         
leg0 = legend('ABK N_X(t)','ABK N_Y(t)','DE N_X(t)','DE N_Y(t)');
set(leg0,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');       

%% Plot State Space: nullclines, direction field, trajectories
cap = agents;                               % ** Set appropriate state-space window **

X_nc = double(subs(X_nc_sym,X,startvalue:cap));      
Y_nc = double(subs(Y_nc_sym,X,startvalue:cap));

[Xm,Ym] = meshgrid(startvalue:40:cap);     % Mesh-grid values for constructing state space

dX_num = double(subs(dX_sym,[X,Y],{Xm,Ym}));
dY_num = double(subs(dY_sym,[X,Y],{Xm,Ym}));
% r = ( dX_num.^2 + dY_num.^2 ).^0.5;
r=1;
%%
figure('Name','State Space','NumberTitle','off');               hold on;
plot(startvalue:cap,X_nc,'b-.',startvalue:cap,Y_nc,'m-.','LineWidth',2);

% quiver(Xm,Ym,dX_num./r,dY_num./r,'g');                  % plot direction field                           
streakarrow(Xm,Ym,dX_num./r,dY_num./r,0.7,1);           % plot direction field

scatter(y_sol(:,1),y_sol(:,2),'+k');       % plot deterministic trajectory
% plot(Sx,Sy,'.r');                           % plot stochastic trajectory

axis([0 cap 0 cap]);
xlabel('X');                    ylabel('Y');                    
legend('X nullcline','Y nullcline','Dir field','DE','ABK');
% legend('Direction field','Deterministic','ABM stochastic');
%% Finish
clear i j q r x y w temp temp1 temp2 temp3 temp4;
toc

%% Notes
% - Works.