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

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;      clc;     tic;                           
rng(1);

global a b c d;

reps = 100;                          % Repeat simulation this many times

maxTime = 2500;                      % Maximum Simulation time (sec)
dt = 0.02;                           % Constant (fixed) time step increment (sec)
t_steps = maxTime / dt;
time = 0:dt:maxTime;

a = 0.0001;                           % 2 X + Y --> 3 X rate constant (3rd order)
b = 0.30;                              % birth rate constant (0th order) for X
c = 0.075;                             % X --> Y conversion rate constant (1st order)
d = 0.02;                             % death rate constant (1st order) for X

% *** Initial population sizes ***
Xi = 0;                        Yi = 0;

agents = 200;                       % Max number of agents allowed

% ** Preallocate memory **
X_all = zeros(reps,t_steps+1);        Y_all = zeros(reps,t_steps+1);         

%% Repeat simulation 'reps' number of times
for n=1:reps
        
    fprintf('.');
    
    % Initialize - Preallocate memory
    Tx = zeros(1,t_steps+1);            Ty = zeros(1,t_steps+1);
%     rxns = zeros(1,t_steps);            
    P_a = zeros(1,t_steps);

    % ******** Initial conditions - Number of X, Y Agents ************
    Tx(1) = Xi;                         Ty(1) = Yi;    
    % ****************************************************************

    tempX = zeros(1,agents);                        
    % Put "1" where agents are "alive", then randomize the array
    for w=1:Tx(1),                      tempX(w)=1;             end
    tempX = RandArray(tempX);           % Randomize array
    Xv = [tempX ; tempX];               % Initialize vector storing state of X agents     

    tempY = zeros(1,agents); 
    % Put "1" where agents are "alive", then randomize the array
    for z=1:Ty(1),                      tempY(z)=1;             end
    tempY = RandArray(tempY);           % Randomize array
    Yv = [tempY ; tempY];               % Initialize vector storing state of B agents
    % Notes on Xv, Yv:
    % - Markov process, so only previous and current time steps needed --> 2 rows:          
    clear w z tempX tempY;

    % *** ABK simulation ***

    % Time-independent probabilities
    P_b = b * dt;                    % Probability of "birth of X" process (0th order)
    P_d = 1 - exp(-d * dt);          % Probability of "death" process (1st order for X)
    P_c = 1 - exp(-c * dt);          % Probability of 1st order rx X --> Y (1st order wrt X)

    for t=1:t_steps

%         P_a(t) = a * Ty(t) * (Tx(t)-1) * dt;          % Prob of 2 X + Y --> 3 X  wrt X
        P_a(t) = a * Tx(t) * (Tx(t)-1) * dt;          % Prob of 2 X + Y --> 3 X  wrt Y
        
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
%             elseif r >= P_d+P_c && r < P_d+P_c+P_a(t)
%                 temp2 = find(Yv(1,:)==1);               % Find Y agents who are "alive"
%                 w = ceil(rand*size(temp2,2));           % Pick Y "alive" agent randomly
%                 Yv(2,temp2(w)) = 0;                     % Y agent dies 
%                 
%                 tempXd = find(Xv(1,:)==0);
%                 s = ceil(rand*size(tempXd,2));           % Pick X "dead" agent randomly
%                 Xv(2,tempXd(s)) = 1;                     % X agent is born
            end
        end
 
        temp3 = find(Yv(1,:)==1);
        for j=1:size(temp3,2)
%             Xalive = Tx(t) - 2*rxns(t);  
%             P_a(t) = a * Xalive * (Xalive-1) * dt;      % Prob of 2 X + Y --> 3 X          
            
            if rand < P_a(t)
                Yv(2,temp3(j)) = 0;                     % Y agent dies
                temp4 = find(Xv(1,:)==0);               % Find X agents who are "dead"
                q = ceil(rand*size(temp4,2));           % Pick X "dead" agent randomly
                Xv(2,temp4(q)) = 1;                     % X agent is born 
%                 rxns(t) = rxns(t) + 1;
            end
        end

        Tx(t+1) = sum(Xv(2,:));                  Ty(t+1) = sum(Yv(2,:));

        Xv(1,:) = Xv(2,:);                       Yv(1,:) = Yv(2,:); 
    
    end         % end "for t" loop
    
    X_all(n,:) = Tx;                Y_all(n,:) = Ty;   

end         % end "for n=1:reps" loop

clear Xv Yv temp* i j;
%% Calculate AVERAGE + SDEV Time Course
avg_X = mean(X_all);                    sdev_X = std(X_all);
avg_Y = mean(Y_all);                    sdev_Y = std(Y_all);

%% Symbolic calculations and numerical solution to DEs
syms X Y positive;

DEchoice = input('Canonical (0) or Agent-Based DEs (1): ');

if DEchoice == 0   
    [t_sol_can , y_sol_can] = ode45(@bruss_can_dif,0:maxTime/1000:maxTime,[Xi ; Yi]);   % Canonical DEs
    % For large oscillations, use 'ode23t' for better integration.
    
    dX_sym = b + a * X^2 * Y - (c + d) * X;
    dY_sym =  c * X - a * X^2 * Y;
    
    X_nc_sym = ((c + d)*X - b) / (a*X^2);       % X nullcline. Given in form Y = expression
    Y_nc_sym = c / (a*X);                       % Y nullcline. Given in form Y = expression

    X_ss_can = b / d;                      Y_ss_can = (c * d) / (b * a);  
    
    disp(['(Canonical) Theoretical steady state X = ' num2str(X_ss_can) ' , Y = ' num2str(Y_ss_can)]);

    % Condition for existence of limit cycle (Canonical DEs):
    crit = d + a * b^2 / d^2;
    if c > crit,    disp('Limit cycle expected (caconical DEs)');      end
    
    startvalue = 1;                             % For making mesh-grid
    
    % Linearize and calculate eigenvalues
    Jac = jacobian([dX_sym,dY_sym],[X,Y]);
    J = eval(subs(Jac,[X,Y],[X_ss_can,Y_ss_can]));          

elseif DEchoice == 1
    [t_sol_ab , y_sol_ab] = ode45(@bruss_ab_dif,0:maxTime/1000:maxTime,[Xi ; Yi]);    % Agent-Based DEs
    
    dX_sym = b + a * X * (X-1) * Y - (c + d) * X;
    dY_sym =  c * X - a * X * (X-1) * Y;
    
    X_nc_sym = ((c + d)*X - b) / (a*X*(X-1));   % X nullcline. Given in form Y = expression
    Y_nc_sym = c / (a*(X-1));                   % Y nullcline. Given in form Y = expression

    X_ss_ab = b / d;                      Y_ss_ab = (c * d) / (a*(b - d));
    
    disp(['(Agent-based) Theoretical steady state X = ' num2str(X_ss_ab) ' , Y = ' num2str(Y_ss_ab)]);

    % Condition for existence of limit cycle (Agent-based DEs)
    crit = ( d + a * b * (b-d) / d^2 ) / ( b / (b-d) );
    if c > crit,    disp('Limit cycle expected (Agent-Based DEs)');    end
    
    startvalue = 2;                             % For making mesh-grid; X ~= 1 in Agent-Based DE
    
    % Linearize and calculate eigenvalues
    Jac = jacobian([dX_sym,dY_sym],[X,Y]);
    J = eval(subs(Jac,[X,Y],[X_ss_ab,Y_ss_ab]));          

else
    disp('Wrong choice!');
end

TrJ = J(1,1) + J(2,2);                  DetJ = det(J);
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

%% Plot [selected] time course
trial = 1;
% ul = max(X_all(trial,:));

fig0 = figure('Name',['Time Course: Trial ' num2str(trial)],...
    'NumberTitle','off','Position',[1 1 500 500]);                          hold on;

% title(['N_{X,i} = ' num2str(Xi) ' , N_{Y,i} = ' num2str(Yi)],...
%     'FontSize',12,'FontName','Times New Roman');

title(['k_a = ' num2str(a) ' , k_b = ' num2str(b) ' , k_c = ' num2str(c) ...
    ' , k_d = ' num2str(d) ' sec^{-1}'], 'FontSize',12,'FontName','Times New Roman');

p_Xs = plot(time(1:1:end),X_all(trial,1:1:end),'+r',...
    'MarkerSize',3,'DisplayName','ABK N_{X}(t)');                                
p_Ys = plot(time(1:10:end),Y_all(trial,1:10:end),'g',...
    'MarkerSize',3,'DisplayName','ABK N_{Y}(t)');

if DEchoice == 0
    p_X_DEcan = plot(t_sol_can,y_sol_can(:,1),'--','Color',[0.75 0 1],...
        'LineWidth',2,'DisplayName','Can DE N_{X}(t)');                                       
    p_Y_DEcan = plot(t_sol_can,y_sol_can(:,2),'--','Color',[0 0.5 0],...
        'LineWidth',2,'DisplayName','Can DE N_{Y}(t)');
    leg0 = legend([p_Xs , p_Ys , p_X_DEcan , p_Y_DEcan]);   % ABK + selected DE trajectories
else
    p_X_DEab = plot(t_sol_ab,y_sol_ab(:,1),'--','Color',[0 0.25 0.75],...
        'LineWidth',2,'DisplayName','Ab DE N_{X}(t)');                                       
    p_Y_DEab = plot(t_sol_ab,y_sol_ab(:,2),'--','Color',[0 0.50 0],...
        'LineWidth',2,'DisplayName','Ab DE N_{Y}(t)');
    leg0 = legend([p_Xs , p_Ys , p_X_DEab , p_Y_DEab]);   % ABK + selected DE trajectories
end

axis tight;                                     
% axis([0 time(end) 0 85]);                               
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');                           hold off;
xlabel('t (sec)');                 ylabel('N(t)'); 

% leg0 = legend([p_Xs , p_Ys]);                       % Only ABK trajectories
% leg0 = legend([p_Xs , p_Ys , p_X_DEcan , p_Y_DEcan , p_X_DEab , p_Y_DEab]);   % ABK + all DE trajectories
set(leg0,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthWest');       

% Create textbox
% annotation(fig0,'textbox',[0.03 0.92 0.05 0.07],'String',{'b)'},'FontSize',12,...
%     'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');

%% Spectral Analysis of 'trial' trajectory
fs = 1/dt;                                          % sampling frequency
% t_osc = 0:dt:totalTime; 
X_osc = X_all(trial,1:end);
Y_osc = Y_all(trial,1:end);          
% plot(t_osc,X_osc,t_osc,Y_osc);     xlabel = ('time (sec)');        ylabel('#');

Osc = [X_osc; Y_osc];
f_th = abs(imag(lambda(1))) / (2*pi);          % Theoretical oscillation frequency (Hz)
period_th = 1 / f_th;                               % Theoretical period (sec)

fig1 = figure('Name',['Spectral Analysis: trial #' num2str(trial)],...
    'NumberTitle','off','Position',[1 1 350 550]);

% do FFT on R, F data
n = 2^(nextpow2(length(Y_osc))+4);                      
f = fs / n * (0:n/2);                   

species = {'\underbar{Species X}', '\underbar{Species Y}'}; 

Z = zeros(size(species,2),n);               Zn = zeros(size(species,2),n);

for q=1:size(species,2)
    Z(q,:) = fft(Osc(q,:),n);
    Zn(q,:) = abs(Z(q,:)/n);
    subplot(size(species,2),1,q);                              hold on;
    title(species{q},'FontName','Times New Roman','FontSize',12,'Interpreter','LaTex');  
    plot(1000*f,Zn(q,1:n/2+1));                                   
    
%     axis([0.5*f_th 2*f_th 0 12*q]);
    xlim([100*f_th 1800*f_th]);
    ylim([0 3*Zn(q,find(abs(f-f_th)<=0.002 & abs(f-f_th)>=0.0005,1))]);


    tempy = [0 50 100];                     tempf = 1000*[f_th f_th f_th];           
    plot(tempf,tempy,'r');                  
    set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
    xlabel('f (mHz)');                      ylabel('Power');
end

% Create textbox
annotation(fig1,'textbox',[0.03 0.92 0.05 0.07],'String',{'b)'},'FontSize',12,...
    'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');

hold off;                          clear q tempy tempf species;

%% Plot AVERAGE Time Course - SINGLE PLOT

fig2 = figure('Name','Avg TC','NumberTitle','off','Position',[1 1 500 500]);        hold on;
% title(['N_{X,i} = ' num2str(Xi) ' , N_{Y,i} = ' num2str(Yi)],...
%     'FontSize',12,'FontName','Times New Roman');                            

p_X = plot(time(1:20:end),avg_X(1:20:end),'r','LineWidth',2,'DisplayName','<N_{X}(t)>_{sim}');   
plot(time(1:20:end),avg_X(1:20:end)+sdev_X(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
plot(time(1:20:end),avg_X(1:20:end)-sdev_X(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);

p_Y = plot(time(1:20:end),avg_Y(1:20:end),'g','LineWidth',2,'DisplayName','<N_{Y}(t)>_{sim}');
% plot(time(1:20:end),avg_Y(1:20:end)+sdev_Y(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
% plot(time(1:20:end),avg_Y(1:20:end)-sdev_Y(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);

if DEchoice == 0
    p_X_DEcan = plot(t_sol_can,y_sol_can(:,1),'--','Color',[0.75 0 1],...
        'LineWidth',2,'DisplayName','Can DE N_X(t)');
    p_Y_DEcan = plot(t_sol_can,y_sol_can(:,2),'--','Color',[0 0.5 0],...
        'LineWidth',2,'DisplayName','Can DE N_{Y}(t)');          
    leg2 = legend([p_X , p_Y , p_X_DEcan , p_Y_DEcan]);
else
    p_X_DEab = plot(t_sol_ab,y_sol_ab(:,1),'--','Color',[0.75 0 0.25],...
        'LineWidth',2,'DisplayName','Ab DE N_{X}(t)');                                       
    p_Y_DEab = plot(t_sol_ab,y_sol_ab(:,2),'--','Color',[0 0.70 0],...
        'LineWidth',2,'DisplayName','Ab DE N_{Y}(t)');
    leg2 = legend([p_X , p_Y , p_X_DEab , p_Y_DEab]);
end

% axis tight;                                         
axis([0 time(end) 0 125]);                                                   
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');                   hold off;
xlabel('t (sec)');            ylabel('N(t)');                                        

set(leg2,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');                                          

% Create textbox
% annotation(fig2,'textbox',[0.03 0.92 0.05 0.07],'String',{'b)'},'FontSize',12,...
%     'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');

%% Plot AVERAGE Time Course - SEPARATE plots

fig1 = figure('Name','Avg TC ; SEPARATE Plots',...
    'NumberTitle','off','Position',[1 1 500 800]);        

subplot(2,1,1);                                                        hold on;         
plot(time(1:20:end),avg_X(1:20:end),'r','LineWidth',3);

if DEchoice == 0
    plot(t_sol_can,y_sol_can(:,1),'--','Color',[0.75 0 1],'LineWidth',2);
else
    plot(t_sol_ab,y_sol_ab(:,1),'--','Color',[0.75 0 1],'LineWidth',2);
end

plot(time(1:20:end),avg_X(1:20:end)+sdev_X(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
plot(time(1:20:end),avg_X(1:20:end)-sdev_X(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
% axis([0 time(end) 0 1.05*N_max]);                                 
% axis tight;  
axis([0 time(end) 0 150]);                       % ** Specify y_max **
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('t (sec)');            ylabel('N_X(t)');                                        
leg2 = legend('<N_X(t)>_{sim}','DE  N_X(t)','Location','NorthEast');
set(leg2,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95]);                                     hold off;

subplot(2,1,2);                                                        hold on;
plot(time(1:20:end),avg_Y(1:20:end),'g','LineWidth',3);

if DEchoice == 0
    plot(t_sol_can,y_sol_can(:,2),'--','Color',[0 0.5 0],'LineWidth',2);
else
    plot(t_sol_ab,y_sol_ab(:,2),'--','Color',[0 0.5 0],'LineWidth',2);  
end

plot(time(1:20:end),avg_Y(1:20:end)+sdev_Y(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
plot(time(1:20:end),avg_Y(1:20:end)-sdev_Y(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
% axis([0 time(end) 0 1.05*N_max]);                                % axis tight; 
axis([0 time(end) 0 150]);                       % ** Specify y_max **
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('t (sec)');            ylabel('N_{Y}(t)');                                        
leg3 = legend('<N_{Y}(t)>_{sim}','DE  N_{Y}(t)','Location','SouthEast');
set(leg3,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95]);                                     hold off;

% Create textboxes
annotation(fig1,'textbox',[0.025 0.900 0.05 0.07],'String',{'a)'},'FontSize',12,...
    'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');
annotation(fig1,'textbox',[0.025 0.425 0.05 0.07],'String',{'b)'},'FontSize',12,...
    'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');

%% Spectral Analysis of the AVG trajectories
fs = 1/dt;                                          % sampling frequency
t_osc = 0:dt:maxTime;                          
X_osc = avg_X(t_osc(1)/dt+1:end);
Y_osc = avg_Y(t_osc(1)/dt+1:end);

% plot(t_osc,X_osc,t_osc,Y_osc);      % xlabel = ('time (sec)');      ylabel('#');

Osc = [X_osc; Y_osc];

figure('Name','Spectral Analysis: AVG','NumberTitle','off','Position',[1 1 350 550]);          

% do FFT on R, Ep, X data
n = 2^(nextpow2(length(Y_osc))+4);                      
f = fs / n * (0:n/2);                   species = {'Species X', 'Species Y'};  

Z  = zeros(size(species,2),n);
Zn = zeros(size(species,2),n);

for q=1:size(species,2)
    Z(q,:) = fft(Osc(q,:),n);
    Zn(q,:) = abs(Z(q,:)/n);
    subplot(size(species,2),1,q);        title(species{q});              hold on;
    plot(1000*f,Zn(q,1:n/2+1));  
    
%     axis([0.5*f_th 2*f_th 0 12*q]);
    xlim([100*f_th 1800*f_th]);
    ylim([0 5*Zn(q,find(abs(f-f_th)<=0.002 & abs(f-f_th)>=0.0005,1))]);


    tempy = [0 50 100];                     tempf = 1000*[f_th f_th f_th];           
    plot(tempf,tempy,'r');                  
    set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
    xlabel('f (mHz)');                      ylabel('Power');
end

hold off;                                   clear q tempy tempf species;

%% Calculate Nullclines, meshgrid for direction field
% cap = 5*agents/10;                           % ** Set appropriate state-space window **
cap = 110;                                   % ** Set appropriate state-space window **

X_nc = double(subs(X_nc_sym,X,startvalue:cap));      
Y_nc = double(subs(Y_nc_sym,X,startvalue:cap));

[Xm,Ym] = meshgrid(startvalue:cap/20:cap);     % Mesh-grid values for constructing state space

dX_num = double(subs(dX_sym,[X,Y],{Xm,Ym}));
dY_num = double(subs(dY_sym,[X,Y],{Xm,Ym}));
% r = ( dX_num.^2 + dY_num.^2 ).^0.5;
r=1;
%% Plot Nullclines, State-Space trajectories
figure('Name','State Space','NumberTitle','off','Position',[1 1 500 500]);          hold on;            

p_ncX = plot(startvalue:cap,X_nc,'Color',[1 0.3 0.2],...
    'LineWidth',2,'DisplayName','Ab N_X nullcline');
p_ncY = plot(startvalue:cap,Y_nc,'Color',[0.4 0.4 1],...
    'LineWidth',2,'DisplayName','Ab N_Y nullcline');            

% quiver(Xm,Ym,dX_num./r,dY_num./r,'g');                  % plot direction field                           
streakarrow(Xm,Ym,dX_num./r,dY_num./r,0.7,1);           % plot direction field

% p_abk = plot(avg_X,avg_Y,'b','LineWidth',2,...
%     'DisplayName','Mean ABK trajectory');           % plot AVG stochastic trajectory
% p_abk = plot(X_all(trial,1:10:end),Y_all(trial,1:10:end),'b','LineWidth',2,...
%     'DisplayName','Trial ABK trajectory');           % plot TRIAL stochastic trajectory

p_ic = plot(Xi,Yi,'rp','MarkerSize',9,'MarkerFaceColor','r',...
    'DisplayName',['Initial: (' num2str(Xi) ' , ' num2str(Yi) ')']);     % plot initial condition

if DEchoice == 0
    p_DEcan = plot(y_sol_can(:,1),y_sol_can(:,2),'Color',[0.7 0.7 0.7],'LineWidth',2,...
        'DisplayName','Can DE trajectory');                 % plot canonical deterministic trajectory
    p_fp_can = plot(X_ss_can,Y_ss_can,'oc','MarkerSize',8,...
        'MarkerFaceColor',[1 0.6 0],'DisplayName',...
        ['Can FP: (' num2str(X_ss_can) ' , ' num2str(Y_ss_can) ')']);  % Canonical FP
    leg = legend([p_ncX, p_ncY, p_fp_can , p_ic, p_DEcan ]);
else
    p_DEab = plot(y_sol_ab(:,1),y_sol_ab(:,2),'k','LineWidth',1,...
        'DisplayName','Ab DE trajectory');                 % plot agent-based deterministic trajectory
    p_fp_ab = plot(X_ss_ab,Y_ss_ab,'oc','MarkerSize',8,...
        'MarkerFaceColor',[0.5 0.8 0.2],'DisplayName',...
        ['Ab FP: (' num2str(X_ss_ab) ' , ' num2str(Y_ss_ab) ')']);    % Agent-based FP
    leg = legend([p_ncX, p_ncY , p_fp_ab , p_ic , p_DEab]);
end

set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
axis([0 cap 0 cap]);            set(gca,'DataAspectRatio',[1 1 1]);            
xlabel('N_X');                  ylabel('N_Y');                  

set(leg,'FontName','Times New Roman','FontSize',8,'Interpreter','TeX',...   
    'EdgeColor',[0.95 0.95 0.95]);                                                  hold off;

%% Finish
% save('matlab.mat','a','b','c','d','agents','Xi','Yi','reps','maxTime','dt','time','*_all');

fprintf(['\n' ElapsedTime(toc)]);

clear R Ep dX* dY* X_nc Y_nc;
clear h i j k m n p q r t v temp*;
clear fig* leg* p_*;
