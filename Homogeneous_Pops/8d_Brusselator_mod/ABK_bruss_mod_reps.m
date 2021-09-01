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

% Or, alternatively stated: 
% 
%   --> X           (birth of X - 0th order; rate constant: b)
% X -->             (death or degradadation of X - 1st order wrt X; rate constant: d)
% X --> Y           (1st order wrt X; rate constant: c)

% Replace the 3rd order reaction
% 2X + Y --> 3X     (3rd order: 2nd wrt X, 1st wrt Y; rate constant: a)
% with: 
% 2X <--> Z                     k_f/k_r = K_eq
% Z + Y --> Z + X               2nd order rate constant a 

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;      clc;     tic;                           
rng(1);

global a b c d k_f k_r;

reps = 100;                          % Repeat simulation this many times

maxTime = 2500;                      % Maximum Simulation time (sec)
dt = 0.01;                           % Constant (fixed) time step increment (sec)
t_steps = maxTime / dt;
time = 0:dt:maxTime;

a = 0.01;                            % Z + Y --> Z + X rate constant (2nd order)
b = 0.3;                             % birth rate constant (0th order) for X
c = 0.075;                           % X --> Y conversion rate constant (1st order)
d = 0.02;                            % death rate constant (1st order) for X
k_f = 0.01;
k_r = 1;

% *** Initial population sizes ***
Xi = 0;                        Yi = 0;                  Zi = 0;        
agents = 200;                       % Max number of agents allowed

% ** Preallocate memory **
X_all = zeros(reps,t_steps+1);      Y_all = zeros(reps,t_steps+1);    Z_all = zeros(reps,t_steps+1);         
%% Repeat simulation 'reps' number of times
for n=1:reps
        
    fprintf('.');
    
    % Initialize - Preallocate memory
    Tx = zeros(1,t_steps+1);            Ty = zeros(1,t_steps+1);        Tz = zeros(1,t_steps+1);
    P_a = zeros(1,t_steps+1);           P_f = zeros(1,t_steps+1); 

    % ******** Initial conditions - Number of X, Y Agents *************************
    Tx(1) = Xi;                         Ty(1) = Yi;                     Tz(1) = Zi;  
    % *****************************************************************************

    tempX = zeros(1,agents);                        
    % Put "1" where agents are "alive", then randomize the array
    for w=1:Tx(1),                      tempX(w)=1;             end
    tempX = RandArray(tempX);           % Randomize array
    Xv = [tempX ; tempX];               % Initialize vector storing state of X agents     

    tempY = zeros(1,agents); 
    % Put "1" where agents are "alive", then randomize the array
    for z=1:Ty(1),                      tempY(z)=1;             end
    tempY = RandArray(tempY);           % Randomize array
    Yv = [tempY ; tempY];               % Initialize vector storing state of Y agents

    tempZ = zeros(1,agents); 
    % Put "1" where agents are "alive", then randomize the array
    for u=1:Tz(1),                      tempZ(u)=1;             end
    tempZ = RandArray(tempZ);           % Randomize array
    Zv = [tempZ ; tempZ];               % Initialize vector storing state of Z agents
    % Notes on Xv, Yv, Zv:
    % - Markov process, so only previous and current time steps needed --> 2 rows:          
    clear w z tempX tempY tempZ;

    % *** ABK simulation ***

    % Time-independent probabilities
    P_b = b * dt;                    % Probability of "birth of X" process (0th order)
    P_d = 1 - exp(-d * dt);          % Probability of "death" process (1st order for X)
    P_c = 1 - exp(-c * dt);          % Probability of 1st order rx X --> Y (1st order wrt X)
    P_r = 1 - exp(-k_r * dt);        % Probability of Z --> 2X (1st order wrt Z)

    for t=1:t_steps
        
        % ****** Time-dependent probabilities  ******
        P_f(t) = k_f * (Tx(t) - 1) * dt;                % Prob of X + X --> Z ,wrt X  (2nd order)
        P_a(t) = a * Tz(t) * dt;                        % Prob of Z + Y --> Z + X ,wrt Y (2nd order)
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
        
        Tx(t+1) = sum(Xv(2,:));           Ty(t+1) = sum(Yv(2,:));       Tz(t+1) = sum(Zv(2,:));

        Xv(1,:) = Xv(2,:);                Yv(1,:) = Yv(2,:);            Zv(1,:) = Zv(2,:); 
    
    end         % end "for t" loop
    
    X_all(n,:) = Tx;                Y_all(n,:) = Ty;                Z_all(n,:) = Tz; 

end         % end "for n=1:reps" loop

clear Xv Yv Zv temp* i j s g h p q;
%% Calculate AVERAGE + SDEV Time Course
avg_X = mean(X_all);                    sdev_X = std(X_all);
avg_Y = mean(Y_all);                    sdev_Y = std(Y_all);
avg_Z = mean(Z_all);                    sdev_Z = std(Z_all);

%% Symbolic calculations and numerical solution to DEs
syms X Y Z positive;
K = k_f / k_r;

DEchoice = input('Canonical (0) or Agent-Based DEs (1):');

if DEchoice == 0   
    [t_sol, y_sol] = ode23t(@bruss_can_dif,0:maxTime/1000:maxTime,[Xi ; Yi ; Zi]);   % Canonical DEs
    
    dX_sym = b - 2 * k_f * X^2 + 2 * k_r * Z + a * Z * Y  - (c + d) * X; 
    dY_sym = c * X - a * Z * Y;
    dZ_sym = k_f * X^2 - k_r * Z;
    
    X_nc_sym = ((c + d)*X - b) / (a*K * X^2);   % X nullcline. Given in form X = expression
    Y_nc_sym = c / (a*K*X);                     % Y nullcline. Given in form Y = expression
    Z_nc_sym = K * X^2;                         % Z nullcline. Given in form Z = expression
    
    X_ss = b / d;                  Y_ss = (c * d) / (a*K * b);         Z_ss = K * b^2 / d^2;
        
    startvalue = 1;                             % For making mesh-grid
    
elseif DEchoice == 1
    [t_sol, y_sol] = ode45(@bruss_ab_dif,0:maxTime/1000:maxTime,[Xi ; Yi ; Zi]);    % Agent-Based DEs
    
    dX_sym = b - 2 * k_f * X * (X-1) + 2 * k_r * Z + a * Z * Y - (c + d) * X;
    dY_sym =  c * X - a * Z * Y;
    dZ_sym = k_f * X * (X-1) - k_r * Z;
    
    X_nc_sym = ((c + d)*X - b) / (a*K * X*(X-1));   % X nullcline. Given in form X = expression
    Y_nc_sym = c / (a*K * (X-1));                   % Y nullcline. Given in form Y = expression
    Z_nc_sym = K * X * (X-1);                       % Z nullcline. Given in form Z = expression
    
    X_ss = b / d;              Y_ss = (c * d) / (a*K * (b - d));         Z_ss = K * b*(b-d) / d^2;
        
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

%% Plot [selected] time course
trial = 1;
% ul = max(X_all(trial,:));

fig0 = figure('Name',['Time Course: Trial ' num2str(trial)],...
    'NumberTitle','off','Position',[1 1 500 500]);                          hold on;

% title(['N_{X,i} = ' num2str(Xi) ' , N_{Y,i} = ' num2str(Yi)],...
%     'FontSize',12,'FontName','Times New Roman');

title(['k_a\prime \cdot K = ' num2str(a*K) ' , k_b = ' num2str(b) ' , k_c = ' num2str(c) ...
    ' , k_d = ' num2str(d) ' sec^{-1}'], 'FontSize',12,'FontName','Times New Roman');

p_Xs = plot(time(1:20:end),X_all(trial,1:20:end),'+r',...
    'MarkerSize',2,'DisplayName','ABK N_{X}(t)');                                
p_Ys = plot(time(1:20:end),Y_all(trial,1:20:end),'g',...
    'MarkerSize',2,'DisplayName','ABK N_{Y}(t)');
p_Zs = plot(time(1:20:end),Z_all(trial,1:20:end),'b',...
    'MarkerSize',2,'DisplayName','ABK N_{Z}(t)');

% p_Xd = plot(t_sol,y_sol(:,1),'--','Color',[0.75 0 1],'LineWidth',2,'DisplayName','DE N_{X}(t)');                                       
p_Yd = plot(t_sol,y_sol(:,2),'--','Color',[0 0.5 0],'LineWidth',2,'DisplayName','DE N_{Y}(t)');
% p_Zd = plot(t_sol,y_sol(:,3),'--','Color',[0 0.5 0.5],'LineWidth',2,'DisplayName','DE N_{Z}(t)');

% axis tight;                                     
axis([0 time(end) 0 125]);                               
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');                           hold off;
xlabel('t (sec)');                 ylabel('N(t)'); 

leg0 = legend([p_Xs , p_Ys , p_Zs]);                          % Only ABK trajectories
% leg0 = legend([p_Xs , p_Ys , p_Zs , p_Xd , p_Yd, p_Zd]);      % ABK + DE trajectories
set(leg0,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');       

% Create textbox
% annotation(fig0,'textbox',[0.03 0.92 0.05 0.07],'String',{'a)'},'FontSize',12,...
%     'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');

%% Spectral Analysis of 'trial' trajectory
fs = 1/dt;                                          % sampling frequency
% t_osc = 0:dt:totalTime; 
X_osc = X_all(trial,1:end);     Y_osc = Y_all(trial,1:end);         Z_osc = Z_all(trial,1:end);   
% plot(t_osc,X_osc,t_osc,Y_osc);     xlabel = ('time (sec)');        ylabel('#');

Osc = [X_osc; Y_osc ; Z_osc];
f_th = abs(imag(lambda(2))) / (2*pi);          % Theoretical oscillation frequency (Hz)
period_th = 1 / f_th;                               % Theoretical period (sec)

fig1 = figure('Name',['Spectral Analysis: trial #' num2str(trial)],...
    'NumberTitle','off','Position',[1 1 350 750]);

% do FFT on R, F data
n = 2^(nextpow2(length(Y_osc))+4);                      
f = fs / n * (0:n/2);                   

species = {'\underbar{Species X}', '\underbar{Species Y}' , '\underbar{Species Z}'}; 

Z = zeros(size(species,2),n);               Zn = zeros(size(species,2),n);

for q=1:size(species,2)
    Z(q,:) = fft(Osc(q,:),n);
    Zn(q,:) = abs(Z(q,:)/n);
    subplot(size(species,2),1,q);                              hold on;
    title(species{q},'FontName','Times New Roman','FontSize',12,'Interpreter','LaTex');  
    plot(1000*f,Zn(q,1:n/2+1));                                   
    
%     axis([0.5*f_th 2*f_th 0 12*q]);
    xlim([100*f_th 1800*f_th]);
    ylim([0 10*Zn(q,find(abs(f-f_th)<=0.002 & abs(f-f_th)>=0.0005,1))]);


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

fig2 = figure('Name','Avg TC','NumberTitle','off','Position',[1 1 500 500]);  
title(['N_{X,i} = ' num2str(Xi) ' , N_{Y,i} = ' num2str(Yi) ' , N_{Z,i} = ' num2str(Zi)],...
    'FontSize',12,'FontName','Times New Roman');                            hold on;

p_X = plot(time(1:20:end),avg_X(1:20:end),'r','LineWidth',3,'DisplayName','<N_{X}(t)>_{sim}');   
plot(time(1:20:end),avg_X(1:20:end)+sdev_X(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
plot(time(1:20:end),avg_X(1:20:end)-sdev_X(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);

p_Y = plot(time(1:20:end),avg_Y(1:20:end),'g','LineWidth',3,'DisplayName','<N_{Y}(t)>_{sim}');
plot(time(1:20:end),avg_Y(1:20:end)+sdev_Y(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
plot(time(1:20:end),avg_Y(1:20:end)-sdev_Y(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);

p_Z = plot(time(1:20:end),avg_Z(1:20:end),'b','LineWidth',3,'DisplayName','<N_{Z}(t)>_{sim}');
plot(time(1:20:end),avg_Z(1:20:end)+sdev_Z(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
plot(time(1:20:end),avg_Z(1:20:end)-sdev_Z(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);

p_Xd = plot(t_sol,y_sol(:,1),'--','Color',[0.75 0 1],'LineWidth',2,'DisplayName','DE N_X(t)');
p_Yd = plot(t_sol,y_sol(:,2),'--','Color',[0 0.5 0],'LineWidth',2,'DisplayName','DE N_{Y}(t)');                                     
p_Zd = plot(t_sol,y_sol(:,3),'--','Color',[0 0.5 0.5],'LineWidth',2,'DisplayName','DE N_{Z}(t)');

% axis tight;                                         
axis([0 time(end) 0 125]);                                                   
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');                   hold off;
xlabel('t (sec)');            ylabel('N(t)');                                        

leg2 = legend([p_X , p_Y , p_Z , p_Xd , p_Yd , p_Zd]);
set(leg2,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95],'Location','NorthEast');                                          

% Create textbox
annotation(fig2,'textbox',[0.03 0.92 0.05 0.07],'String',{'b)'},'FontSize',12,...
    'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');

%% Plot AVERAGE Time Course - SEPARATE plots

fig1 = figure('Name','Avg TC ; SEPARATE Plots',...
    'NumberTitle','off','Position',[1 1 500 1200]);        

subplot(3,1,1);                                                        hold on;         
plot(time(1:20:end),avg_X(1:20:end),'r','LineWidth',3);
plot(t_sol,y_sol(:,1),'--','Color',[0.75 0 1],'LineWidth',2);   
plot(time(1:20:end),avg_X(1:20:end)+sdev_X(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
plot(time(1:20:end),avg_X(1:20:end)-sdev_X(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
% axis([0 time(end) 0 1.05*N_max]);                                 
% axis tight;  
axis([0 time(end) 0 70]);                       % ** Specify y_max **
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('t (sec)');            ylabel('N_X(t)');                                        
leg2 = legend('<N_X(t)>_{sim}','DE  N_X(t)','Location','NorthEast');
set(leg2,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95]);                                     hold off;

subplot(3,1,2);                                                        hold on;
plot(time(1:20:end),avg_Y(1:20:end),'g','LineWidth',3);
plot(t_sol,y_sol(:,2),'--','Color',[0 0.5 0],'LineWidth',2);  
plot(time(1:20:end),avg_Y(1:20:end)+sdev_Y(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
plot(time(1:20:end),avg_Y(1:20:end)-sdev_Y(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
% axis([0 time(end) 0 1.05*N_max]);                                % axis tight; 
axis([0 time(end) 0 125]);                       % ** Specify y_max **
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('t (sec)');            ylabel('N_{Y}(t)');                                        
leg3 = legend('<N_{Y}(t)>_{sim}','DE  N_{Y}(t)','Location','NorthEast');
set(leg3,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95]);                                     hold off;

subplot(3,1,3);                                                        hold on;
plot(time(1:20:end),avg_Z(1:20:end),'b','LineWidth',3);
plot(t_sol,y_sol(:,3),'--','Color',[0 0.5 0.5],'LineWidth',2);  
plot(time(1:20:end),avg_Z(1:20:end)+sdev_Z(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
plot(time(1:20:end),avg_Z(1:20:end)-sdev_Z(1:20:end),'LineStyle','--','Color',[0.8 0.8 0.8]);
% axis([0 time(end) 0 1.05*N_max]);                                % axis tight; 
axis([0 time(end) 0 25]);                       % ** Specify y_max **
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('t (sec)');            ylabel('N_{Z}(t)');                                        
leg3 = legend('<N_{Z}(t)>_{sim}','DE  N_{Z}(t)','Location','NorthEast');
set(leg3,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95]);                                     hold off;

% Create textboxes
% annotation(fig1,'textbox',[0.025 0.900 0.05 0.07],'String',{'a)'},'FontSize',12,...
%     'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');
% annotation(fig1,'textbox',[0.025 0.425 0.05 0.07],'String',{'b)'},'FontSize',12,...
%     'FontName','Times New Roman','FitBoxToText','off','LineStyle','none');

%% Spectral Analysis of the AVG trajectories
fs = 1/dt;                                          % sampling frequency
t_osc = 0:dt:maxTime;                          
X_osc = avg_X(t_osc(1)/dt+1:end);     Y_osc = avg_Y(t_osc(1)/dt+1:end);     
Z_osc = avg_Z(t_osc(1)/dt+1:end);

% plot(t_osc,X_osc,t_osc,Y_osc);      % xlabel = ('time (sec)');      ylabel('#');

Osc = [X_osc; Y_osc; Z_osc];

figure('Name','Spectral Analysis: AVG','NumberTitle','off','Position',[1 1 350 750]);          

% do FFT on R, Ep, X data
n = 2^(nextpow2(length(Y_osc))+4);                      
f = fs / n * (0:n/2);                   species = {'Species X', 'Species Y', 'Species Z'};  

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

%% Calculate meshgrid for direction field ; Plot State-Space trajectories
% cap = max(max(y_sol))*1.2;                          % ** Set appropriate state-space window **

% [Xm,Ym,Zm] = meshgrid(startvalue:cap/10:cap);       % Mesh-grid values for constructing state space

% dX_num = double(subs(dX_sym,[X,Y,Z],{Xm,Ym,Zm}));
% dY_num = double(subs(dY_sym,[X,Y,Z],{Xm,Ym,Zm}));
% dZ_num = double(subs(dZ_sym,[X,Y,Z],{Xm,Ym,Zm}));

% r = ( dX_num.^2 + dY_num.^2 + dZ_num.^2 ).^0.5;
% r=1;

% figure('Name','State Space','NumberTitle','off','Position',[1 1 500 500]);          hold on;            
% 
% quiver3(Xm,Ym,Zm,dX_num./r,dY_num./r,dZ_num./r,'g');      % plot direction field                           
% 
% p_de = plot3(y_sol(:,1),y_sol(:,2),y_sol(:,3),'k','LineWidth',1,...
%     'DisplayName','DE trajectory');                 % plot deterministic trajectory
% p_abk = plot3(avg_X,avg_Y,avg_Z,'b','LineWidth',2,...
%     'DisplayName','Mean ABK trajectory');           % plot AVG stochastic trajectory
% p_abk = plot3(X_all(trial,1:10:end),Y_all(trial,1:10:end),Z_all(trial,1:10:end),'b',...
%     'LineWidth',2,'DisplayName','Trial ABK trajectory');       % plot TRIAL stochastic trajectory

% p_ic = plot3(Xi,Yi,Zi,'rp','MarkerSize',9,'MarkerFaceColor','r',...  % plot initial condition
%     'DisplayName',['Initial: (' num2str(Xi) ',' num2str(Yi) ',' num2str(Zi) ')']);
% p_fp = plot3(X_ss,Y_ss,Z_ss,'oc','MarkerSize',8,...
%     'MarkerFaceColor','y','DisplayName',...
%     ['FP: (' num2str(X_ss) ',' num2str(Y_ss) ',' num2str(Z_ss) ')']);             % FP
% 
% set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
% axis([0 cap 0 cap 0 cap]);            set(gca,'DataAspectRatio',[1 1 1]);            
% xlabel('N_X');                  ylabel('N_Y');                  zlabel('N_Z');
% 
% leg = legend([p_abk , p_de , p_fp , p_ic]);
% set(leg,'FontName','Times New Roman','FontSize',8,'Interpreter','TeX',...   
%     'EdgeColor',[0.95 0.95 0.95]);                                                  hold off;

%% Finish
% save('matlab.mat','a','b','c','d','k_*','agents','*i','reps','maxTime','dt','time','*_all');

fprintf(['\n' ElapsedTime(toc)]);

clear R Ep dX* dY* X_nc Y_nc;
clear h i j k m n p q r t v temp*;
clear fig* leg* p_*;
