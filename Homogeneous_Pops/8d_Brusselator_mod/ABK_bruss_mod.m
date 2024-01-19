% The Brusselator (a hypothetical system) reaction scheme (reference?). 
% ABK implementation.

%        ^                               b  = birth (0th order) rate constant for X
%      d |                               d  = death (1st order) rate constant for X
%        |       k_f                     
%   -->  X + X  <--->  Z                 k_f/k_r: rate constants for 2X <--> Z
%    b   ^       k_r   |                     
%        |a            |                 a =  Z + Y --> Z + X rate constant (2nd order)
%        | ------------- 
%    c  \/                               c  = X --> Y (1st order) rate constant 
%        Y

% Or, alternately stated: 
% 
%   --> X           (birth of X - 0th order; rate constant: b)
% X -->             (death or degradadation of X - 1st order wrt X; rate constant: d)
% X --> Y           (1st order wrt X; rate constant: c)

% Replace the 3rd order reaction
% 2X + Y --> 3X     (3rd order: 2nd wrt X, 1st wrt Y; rate constant: a)
% with: 
% 2X <--> Z                     k_f/k_r = K_eq
% Z + Y --> Z + X               2nd order rate constant a 

% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear;          clc;         tic;
rng(0);
% Declare variables and functions
global a b c d k_f k_r;

maxTime = 2000;                          % Maximum Simulation time (sec)
dt = 0.01;                              % Constant (fixed) time step increment (sec)
t_steps = maxTime / dt;

a = 0.0005;                               % Z + Y --> Z + X rate constant (2nd order)
b = 0.3;                               % birth rate constant (0th order) for X
c = 0.075;                                % X --> Y conversion rate constant (1st order)
d = 0.05;                               % death rate constant (1st order) for X
k_f = 0.1;
k_r = 0.1;
K = k_f / k_r;

Xi = 0;         Yi = 0;         Zi = 0;         % ** Initial condition **

agents = 200;                       % disp(['Agents = ' num2str(agents)]);

%% Initialize - Preallocate memory for variables; ** Initial conditions **
Sx = zeros(1,t_steps);              % Sum of X agents in each time step	
Sy = zeros(1,t_steps);              % Sum of Y agents in each time step
Sz = zeros(1,t_steps);              % Sum of Z agents in each time step

% ******** Initial conditions - Number of A, B agents ********
Sx(1) = Xi;               Sy(1) = Yi;           Sz(1) = Zi;					
% ************************************************************

tempX = zeros(1,agents/2);                        
% Put "1" where agents are "alive", then randomize the array
for w=1:Sx(1),                      tempX(w)=1;             end
tempX = RandArray(tempX);           % Randomize array
Xv = [tempX ; tempX];               % Initialize vector storing state of X agents     

tempY = zeros(1,agents); 
% Put "1" where agents are "alive", then randomize the array
for z=1:Sy(1),                      tempY(z)=1;             end
tempY = RandArray(tempY);           % Randomize array
Yv = [tempY ; tempY];               % Initialize vector storing state of Y agents

tempZ = zeros(1,agents/4); 
% Put "1" where agents are "alive", then randomize the array
for u=1:Sz(1),                      tempZ(u)=1;             end
tempZ = RandArray(tempZ);           % Randomize array
Zv = [tempZ ; tempZ];               % Initialize vector storing state of Z agents

% Notes on Xv, Yv, Zv:
% - Markov process, so only previous and current time steps needed --> 2 rows:          
clear w z u tempX tempY tempZ;

%% ABK Simulation
P_a = zeros(1,t_steps);          % Preallocate memory
P_f = zeros(1,t_steps);          % Preallocate memory

% Time-independent probabilities
P_b = b * dt;                    % Probability of "birth of X" process (0th order)
P_d = 1 - exp(-d * dt);          % Probability of "death" process (1st order for X)
P_c = 1 - exp(-c * dt);          % Probability of 1st order rx X --> Y (1st order wrt X)
P_r = 1 - exp(-k_r * dt);        % Probability of Z --> 2X (1st order wrt Z)

for t=1:t_steps
    
    % ****** Time-dependent probabilities  ******
    P_f(t) = k_f * (Sx(t) - 1) * dt;                % Prob of X + X --> Z ,wrt X  (2nd order)
    P_a(t) = a * Sz(t) * dt;                        % Prob of Z + Y --> Z + X ,wrt Y (2nd order)
    % ****** ****************************  ******

    tempXa = find(Xv(1,:)==1);                      % find "alive" X agents
    tempXd = find(Xv(1,:)==0);                      % find "dead" X agents
    tempYa = find(Yv(1,:)==1);                      % find Y agents who are in "alive" state
    tempYd = find(Yv(1,:)==0);                      % find "dead" Y agents
    tempZa = find(Zv(1,:)==1);                      % find Z agents who are in "alive" state
    tempZd = find(Zv(1,:)==0);                      % Find Z agents who are in "dead" state

    % Take care of 0th order process first    
    if rand < P_b                                   % "Birth", 0th order reaction        
        q = ceil(rand * size(tempXd,2));            % pick a dead X agent randomly
        Xv(2,tempXd(q)) = 1;                        % an X agent is born
        tempXd(q) = 0;      % making sure this X agent is not 'reborn' in this time step
                            % through rx Z --> X + X and Z + Y --> Z + X
    end                                    
    
    % Degradation of X    or    X --> Y     or     2X --> Z    
    for i=1:size(tempXa,2)
        % Skip agents which "die" in this time step due to the rx 2X --> Z
        if tempXa(i) == 0,       continue,           end   
        r = rand;
        if r < P_d
            
            Xv(2,tempXa(i)) = 0;                    % X agent dies
            
        elseif r >= P_d && r < P_d + P_c
            Xv(2,tempXa(i)) = 0;                    % X agent dies
            tempXa(i) = 0;
            
            p = ceil(rand*size(tempYd,2));          % Pick Y "dead" agent randomly
            Yv(2,tempYd(p)) = 1;                    % Y agent is born
            
        elseif r >= P_d + P_c && r < P_d + P_c + P_f(t)
                        
            f = ceil(rand*size(tempXa,2));           % Pick *another* X "alive" agent randomly
            count = 0;
            while f==i || tempXa(f) == 0             % make sure this X agent hasn't already ...
                f = ceil(rand*size(tempXa,2));       % died in this time step.
                count = count + 1;
                if count == 100,        break;      end     % break the 'while' loop
            end  
            if count == 100,        continue;      end     % skip this turn in the 'for' loop
            Xv(2,tempXa(i)) = 0;                     % 1st X agent dies
            tempXa(i) = 0; 
            Xv(2,tempXa(f)) = 0;                     % 2nd X agent dies            
            tempXa(f) = 0;   % making sure this X agent is not sampled again as "alive" in this time step
            
            k = ceil(rand*size(tempZd,2));           % Pick Z "dead" agent randomly
            Zv(2,tempZd(k)) = 1;                     % Z agent is born
            
        end
    end
        
    for h=1:size(tempZa,2)
        if rand < P_r                               % Z --> 2X (1st order wrt Z)
            Zv(2,tempZa(h)) = 0;                    % Z agent dies
                        
            g = ceil(rand(1,2)*size(tempXd,2));     % Pick two X "dead" agents randomly
            while g(1) == g(2) || isempty(find(tempXd(g)==0, 1))==0  
                % Make sure two X agents are different/distinct
                % and haven't already been "born" during
                % this time step in 0th order rx (--> X)
                % (or 2nd order Z + Y --> Z + X, in principle)
                g = ceil(rand(1,2)*size(tempXd,2));
            end
            Xv(2,tempXd(g(1))) = 1;                 % 1st X agent is born
            Xv(2,tempXd(g(2))) = 1;                 % 2nd X agent is born
            
            tempXd(g) = 0;      % making sure these two X agents are not 'reborn' in this 
                                % time step through rxs:  null --> X and Z + Y --> Z + X
        end
    end
        
    for j=1:size(tempYa,2)
        if rand < P_a(t)
            Yv(2,tempYa(j)) = 0;                    % Y agent dies
            
            s = ceil(rand*size(tempXd,2));          % Pick X "dead" agent randomly
            while tempXd(s) == 0                    % make sure this X agent hasn't already
                s = ceil(rand*size(tempXd,2));      % been born in this time step.
            end
            Xv(2,tempXd(s)) = 1;                    % X agent is born 
            
            tempXd(s) = 0;      % making sure this X agent is not 'reborn' in this 
                                % time step through rx  null --> X and Z --> X + X
            % This last declaration is unnecessary since it is the end of the 
            % algorithm in this time step. Still, I am including it for completeness.
        end
    end
            
    Sx(t+1) = sum(Xv(2,:));        Sy(t+1) = sum(Yv(2,:));        Sz(t+1) = sum(Zv(2,:));
    
    Xv(1,:) = Xv(2,:);             Yv(1,:) = Yv(2,:);             Zv(1,:) = Zv(2,:); 
    
end

%% Symbolic calculations and numerical solution to DEs
syms X Y Z positive;

DEchoice = input('Canonical (0) or Agent-Based DEs (1):');

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

%% Graph ABK (stochastic) and deterministic results
figure('Name','Sample Time course','NumberTitle','off');
time = 0:dt:maxTime;                                hold on;
plot(time,Sx,'b');                                        
plot(time,Sy,'r');
plot(time,Sz,'g');
plot(t_sol,y_sol(:,1),'--c','LineWidth',3);                  
plot(t_sol,y_sol(:,2),'--m','LineWidth',3);                               
plot(t_sol,y_sol(:,3),'--y','LineWidth',3);                               

% axis([0 t_sol(end) 0 agents]);          
axis tight;
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');           hold off;

xlabel('time');         ylabel('# Agents');         
leg0 = legend('ABK N_X(t)','ABK N_Y(t)','ABK N_Z(t)','DE N_X(t)','DE N_Y(t)','DE N_Z(t)');
set(leg0,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');       

%% Calculate meshgrid for direction field ; Plot State-Space trajectories
cap = max(max(y_sol))*1.2;                          % ** Set appropriate state-space window **

[Xm,Ym,Zm] = meshgrid(startvalue:cap/10:cap);       % Mesh-grid values for constructing state space

dX_num = double(subs(dX_sym,[X,Y,Z],{Xm,Ym,Zm}));
dY_num = double(subs(dY_sym,[X,Y,Z],{Xm,Ym,Zm}));
dZ_num = double(subs(dZ_sym,[X,Y,Z],{Xm,Ym,Zm}));

% r = ( dX_num.^2 + dY_num.^2 + dZ_num.^2 ).^0.5;
r=1;

figure('Name','State Space','NumberTitle','off','Position',[1 1 500 500]);          hold on;            

quiver3(Xm,Ym,Zm,dX_num./r,dY_num./r,dZ_num./r,'g');      % plot direction field                           
% streakarrow(Xm,Ym,dX_num./r,dY_num./r,0.7,1);           % plot direction field (2D only)

p_abk = plot3(Sx,Sy,Sz,'.r');                           % plot stochastic trajectory
p_de = plot3(y_sol(:,1),y_sol(:,2),y_sol(:,3),'k','LineWidth',1,...
    'DisplayName','DE trajectory');                      % plot deterministic trajectory

p_ic = plot3(Xi,Yi,Zi,'rp','MarkerSize',9,'MarkerFaceColor','r',...  % plot initial condition
    'DisplayName',['Initial: (' num2str(Xi) ',' num2str(Yi) ',' num2str(Zi) ')']);
p_fp = plot(X_ss,Y_ss,'oc','MarkerSize',8,...
    'MarkerFaceColor','y','DisplayName',...
    ['FP: (' num2str(X_ss) ',' num2str(Y_ss) ',' num2str(Z_ss) ')']);             % FP

set(gca,'XMinorTick','on','YMinorTick','on','ZMinorTick','on','Box','off');
axis([0 cap 0 cap 0 cap]);            set(gca,'DataAspectRatio',[1 1 1]);            
xlabel('N_X');                  ylabel('N_Y');                  zlabel('N_Z');

leg = legend([p_abk , p_de , p_fp , p_ic]);
set(leg,'FontName','Times New Roman','FontSize',8,'Interpreter','TeX',...   
    'EdgeColor',[0.95 0.95 0.95]);                                                  hold off;
%% Finish
clear i j q r x y w temp temp1 temp2 temp3 temp4;
toc
