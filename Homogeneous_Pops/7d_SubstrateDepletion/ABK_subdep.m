%          Ep  <--->  E          Rates: k_f, k_r (Reverse rx activated by R)
%          |     ^               Michaelis-Menten constants: Km_f, Km_r
%           \     \              k_s: synthesis of X
%    k_s     \     |             k_b1,2: synthesis of R
%  S ---> X  --->  R  ------>    k_d: degradation of R                
%           k_b1,2      k_d               ** S: Signal, R: Response **  

% Simulating substrate depletion oscillator process: 
% k_s: X synthesis, 0th order process, UPregulated by S
% k_b1: X --> R, 1st order process wrt X
% k_b2: X --> R, 2nd order process wrt X, Ep (Ep is not consummed)
% k_d: R degradation, 1st order process wrt R
% k_f: E synthesis, MM process
% k_r: Ep synthesis, MM process wrt R (but R is not consummed)
% Km_f, Km_r: MM constants for forward and reverse rxs

% Note: Assume number of S molecules/agents is NOT changing (eg, S is DNA molecules)

% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear;      clc;     tic;
rng(0);
global agents k_s k_b1 k_b2 k_d k_f k_r Km_f Km_r S;

totalTime = 5000;             % Simulation time (sec)
dt = 1/50;                    % Constant (fixed) time step increment (sec)
t_steps = totalTime / dt;
% time = zeros(1,t_steps);

agents = 100;
k_s = 0.075;                         % X synthesis, 0th order process, UPregulated by S
k_b1 = 0.001;                       % X --> R, 1st order process wrt X
k_b2 = 0.003;                       % X --> R, 2nd order process wrt X, Ep (Ep is not consummed)
k_d = 0.04;                         % R degradation, 1st order process wrt R
k_f = 0.5;                        % Ep --> E forward rate
k_r = 0.0125;                        % Ep <-- E reverse rate, UPregulated by R
Km_f = 5;                      % MICROSCOPIC Michaelis-Menten constant for forward rx
Km_r = 5;                      % MICROSCOPIC Michaelis-Menten constant for reverse rx
S = 10;                            % Assume number of S molecules/agents is NOT changing

% Initial population sizes
Ri = 0;                 Epi = 0;                Xi = 0;

Tr = zeros(1,t_steps);              Tx = zeros(1,t_steps);
Tep = zeros(1,t_steps);             Te = zeros(1,t_steps);

P_f = zeros(1,t_steps);             P_r = zeros(1,t_steps);
P_b2 = zeros(1,t_steps);            

% ******** Initial conditions - Number of R, Ep/E, X Agents ********
Tr(1) = Ri;                       Tx(1) = Xi;    
Tep(1) = Epi;                     Te(1) = agents - Epi;					
% ****************************************************************

tempR = zeros(1,10*agents);        tempEp = zeros(1,agents); 
tempX = zeros(1,15*agents);
% Put "1" where agents are "alive", then randomize the array
for c=1:Tr(1),                  tempR(c)=1;             end
for d=1:Tep(1),                 tempEp(d)=1;            end
for e=1:Tx(1),                  tempX(e)=1;             end
tempR = RandArray(tempR);                   % Randomize R array
tempEp = RandArray(tempEp);                 % Randomize Ep array
tempX = RandArray(tempX);                   % Randomize X array
% Markov process, so only previous and current time steps needed --> 2 rows:    
R = [tempR ; tempR];                        % Initialize vector storing state of R agents
Ep = [tempEp ; tempEp];                     % Initialize vector storing state of Ep agents
E = ~ Ep;                                   % Initialize vector storing state of E agents
% - Ep and E are complementary (E + Ep = agents)
X = [tempX ; tempX];                        % Initialize vector storing state of X agents
clear c d e tempR tempEp tempX;

%% ABK simulation
t = 1;

P_s = k_s * S * dt;             % Probability of X synthesis process (0th order)
P_b1 = k_b1 * dt;               % Probability X-->R (1st order wrt X)
P_d = k_d * dt;                 % Probability of R degradation process (1st order wrt R)

while t <= t_steps % && Tep(t)<agents
    
    P_b2(t) = k_b2 * Tep(t) * dt;   % Probability (wrt X) of R synthesis (2nd order)
    P_f(t) = k_f / (Km_f + Tep(t)) * dt;          % wrt each Ep molecule
    P_r(t) = k_r * Tr(t) / (Km_r + Te(t)) * dt;   % wrt each E molecule
    
    % Take care of 0th order processes first    
    if rand < P_s
        tempX = find(X(1,:)==0);    % Randomly choose X agent synthesis 
        X(2,tempX(ceil(rand * size(tempX,2)))) = 1;
    end
    % End of 0th order processes
        
    tempR = find(R(1,:)==1);
    for i = 1:size(tempR,2)
        if rand < P_d                   % Degradation of R
            R(2,tempR(i)) = 0;          % R agent is degraded
        end
    end

    tempE = find(E(1,:)==1);
    for j=1:size(tempE,2)
        if rand < P_r(t)
            E(2,tempE(j)) = 0;          % Conversion of E to Ep
            tempEp = find(Ep(1,:)==0);  % Randomly choose Ep agent synthesis
            Ep(2,tempEp(ceil(rand * size(tempEp,2)))) = 1;
        end
    end
    
    tempEp = find(Ep(1,:)==1);
    for k=1:size(tempEp,2)
        if rand < P_f(t)
            Ep(2,tempEp(k)) = 0;        % Conversion of Ep to E
            tempE = find(E(1,:)==0);    % Randomly choose E agent synthesis
            E(2,tempE(ceil(rand * size(tempE,2)))) = 1;
        end
    end
    
    P_b = P_b1 + P_b2(t);               % Total Probability of X --> R
    tempX = find(X(1,:)==1);
    for m = 1:size(tempX,2)
        if rand < P_b                   % X-->R, from 2 separate processes
            X(2,tempX(m)) = 0;          % X agent is converted to R
            tempR = find(R(1,:)==0);    % Randomly choose R agent synthesis 
            R(2,tempR(ceil(rand * size(tempR,2)))) = 1;            
        end
    end  
    
    Tr(t+1) = sum(R(2,:));              Tx(t+1) = sum(X(2,:));          
    Tep(t+1) = sum(Ep(2,:));            Te(t+1) = sum(E(2,:));
    
    R(1,:) = R(2,:);                    X(1,:) = X(2,:);
    Ep(1,:) = Ep(2,:);                  E(1,:) = E(2,:);
    
    t = t + 1;
end

% Remove unnecessary terminal 0's from arrays
if t < t_steps
    Tr = Tr(1:t);                               Tx = Tx(1:t);                             
    Tep = Tep(1:t);                             Te = Te(1:t);             
    P_f = P_f(1:t);          P_r = P_r(1:t);    P_b2 = P_b2(1:t);
end
clear R Ep E X;
%% Solve DE
if exist('t','var')==0,    finaltime = totalTime;      
else                       finaltime = (t-1) * dt;         end
[t_sol, y_sol] = ode45(@subdep2_dif,0:finaltime/1000:finaltime,[Ri ; Epi ; Xi]);
%% Plot time course
time = 0:dt:finaltime;
figure('Name','Time Course','NumberTitle','off');                           hold on;

plot(time(1:20:end),Tr(1:20:end),'r');                                
plot(time(1:20:end),Tep(1:20:end),'g');
plot(time(1:20:end),Tx(1:20:end),'b');
plot(t_sol,y_sol(:,1),'--','Color',[0.75 0 1],'LineWidth',2);                                       
plot(t_sol,y_sol(:,2),'--','Color',[0 0.5 0],'LineWidth',2);
plot(t_sol,y_sol(:,3),'--','Color',[0.15 1 0.75],'LineWidth',2);                                     

% axis([0 finaltime 0 200]);
axis([0 finaltime 0 max(Tx)]);                               
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');                   hold off;
xlabel('t (sec)');                 ylabel('N(t)');      
leg = legend('ABK N_R(t)','ABK N_{Ep}(t)','ABK N_X(t)','DE N_R(t)','DE N_{Ep}(t)','DE N_X(t)');
set(leg,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95],'Location','Best');                                              

%% Symbolic calculations
syms R Ep X positive;
dR_sym = + k_b1 * X + k_b2 * Ep * X - k_d * R;
dEp_sym = + k_r * (agents - Ep) * R / (Km_r + (agents - Ep)) - k_f * Ep / (Km_f + Ep);
dX_sym = + k_s * S - k_b1 * X - k_b2 * Ep * X;

% ** Numerical solution of steady-state values **
SS = vpasolve([dR_sym,dEp_sym,dX_sym],[R,Ep,X]);    % Num solve foe equns = 0
R_ss = double(SS.R);
Ep_ss = double(SS.Ep);
X_ss = double(SS.X);
% ***********************************************

% ** Symbolic calculation of steady-state values **
% ** Nullclines **
% R_nc_sym = solve(dR_sym == 0,R);         % depends on Ep, X
% Ep_nc_sym = solve(dEp_sym == 0,Ep);      % 2 entries: +/-, depends on R
% X_nc_sym = solve(dX_sym == 0,X);         % depends on Ep
% 
% tempX = subs(X_nc_sym,Ep,Ep_nc_sym(1));
% tempR = subs(R_nc_sym,[Ep,X],[Ep_nc_sym(1),tempX]);
% R_ss = double(solve(tempR==R,R));
% Ep_ss = double(subs(Ep_nc_sym(1),R,R_ss));
% X_ss = double(subs(X_nc_sym,Ep,Ep_ss));
% ***********************************************

fprintf(1,'Steady-States:\t\t\tEigenvalues\n');
Jac = jacobian([dR_sym,dEp_sym,dX_sym],[R,Ep,X]);

for w=1:size(R_ss,1)
    J = double(subs(Jac,[R,Ep,X],[R_ss(w),Ep_ss(w),X_ss(w)]));
    lambdas = eig(J);
    fprintf(['R=' num2str(R_ss(w)) ', Ep=' num2str(Ep_ss(w)) ', X=' num2str(X_ss(w)) ':\n']);
    disp([num2str(lambdas)]);
end
clear w;
%% Plot Nullclines, State-space trajectories
% array = 0:agents; 
% exclude = uint16(k_f/k_r);
% if mod(exclude,1) == 0
%     if isempty(find(array==exclude)) == 0
%         array(exclude + 1) = [];   % Remove R = k_f/k_r entry (avoid div by 0 when calc Ep_nc)
%     end
% end
% 
% R_ofX_nc = subs(R_nc_sym,Ep,Ep_nc_sym(1));
% R_ofX_sym = solve(R_ofX_nc==R,R);
% R_nc  = double(subs(R_ofX_sym(3),X,{array}));
% X_nc = double(subs(X_nc_sym,R,array));

% 
% if size(R_nc,1)>1                % Nullcline values must be positive!
%     for j=size(R_nc,1):-1:1
%         if isempty(find(sign(R_nc(j,:)) == -1, 1,'first')) == false
%             R_nc(j,:) = [];
%         end
%     end
% end
% 
% if size(Ep_nc,1)>1               % Nullcline values must be positive!
%     for k=size(Ep_nc,1):-1:1
%         if isempty(find(sign(Ep_nc(k,:)) == -1, 1,'first')) == false
%             Ep_nc(k,:) = [];
%         end
%     end
% end
% 
% [X_int R_int] = intersections(X_nc,array,array,R_nc);
%% State-space 
% [Rm,Epm] = meshgrid(0:agents/10:agents);     % Mesh-grid values for constructing state space
% Em = agents - Epm;
% dR  = + k_b .* Epm + k_s * S - k_d .* Rm;
% dEp = + k_r .* Em .* Rm ./ (Km_r + Em) - k_f .* Epm ./ (Km_f + Epm);
% r = ( dR.^2 + dEp.^2 ).^0.5;
% % r=1;
% 
figure('Name','State Space','NumberTitle','off');               hold on;
% plot(X_nc,array,'b-.',array,R_nc,'m-.','LineWidth',1);
% quiver(Rm,Epm,dR./r,dEp./r,'g');                                    
scatter(y_sol(:,3),y_sol(:,1),'+k');       % plot deterministic trajectory
scatter(Tx,Tr,3,'.r');                    % plot stochastic trajectory
% plot(X_int,R_int,'ob','MarkerSize',6);
axis([0 2*agents 0 1*agents]);
xlabel('X');                    ylabel('R');                  hold off; 
% legend('R nullcline','Ep nullcline','Direction field','Deterministic','ABM stochastic');
%% Spectral Analysis
fs = 1/dt;                                          % sampling frequency
t_osc = 500:dt:totalTime;                          X_osc = Tx(t_osc(1)/dt+1:end);
Ep_osc = Tep(t_osc(1)/dt+1:end);                   R_osc = Tr(t_osc(1)/dt+1:end);
% figure;     plot(t_osc,X_osc,t_osc,Yp_osc,t_osc,Rp_osc);
% xlabel = ('time (sec)');        ylabel('#');
Osc = [R_osc; Ep_osc; X_osc];
f_th = abs(imag(lambdas(2))) / (2*pi);              % Theoretical oscillation frequency (Hz)
period_th = 1 / f_th;                               % Theoretical period (sec)

figure('Position',[1 1 350 800]);                   
% do FFT on R, Ep, X data
n = 2^(nextpow2(length(X_osc))+4);                      
f = fs / n * (0:n/2);                   species = {'$N_R$', '$N_{Ep}$', '$N_X$'}; 

Z = zeros(size(species,2),n);
Zn = zeros(size(species,2),n);

for b=1:size(species,2)
    Z(b,:) = fft(Osc(b,:),n);
    Zn(b,:) = abs(Z(b,:)/n);
    subplot(3,1,b);                                             hold on;
    title(species{b},'FontName','Times New Roman','FontSize',12,'Interpreter','LaTex');  
    plot(1000*f,Zn(b,1:n/2+1));                                   
    
%     axis([0.5*f_th 2*f_th 0 12*b]);
%     xlim([100*f_th 1800*f_th]);
    xlim([0 6]);                % Adjust this when Im(lambdas)=0 
    ylim([0 1.5*Zn(b,find(abs(f-f_th)<=0.0002 & abs(f-f_th)>=0.00005,1))]);


    tempy = [0 50 100];                     tempf = 1000*[f_th f_th f_th];           
    plot(tempf,tempy,'r');                  
    set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
    xlabel('f (mHz)');                      ylabel('Power');
end
  
hold off;                               
% clear b tempy tempf species;

%% Finish
clear array R Ep dR dEp R_nc Ep_nc R_nc_sym Ep_nc_sym;
clear h i j k n r t temp tempR tempEp tempE;
toc
