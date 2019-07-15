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

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;          % clc;
global agents k_s k_b1 k_b2 k_d k_f k_r Km_f Km_r S;

totalTime = 15000;
dt = 0.1;

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

%% Solve DE
% [t_sol, y_sol] = ode45(@subdep_dif,0:500,[0 ; 0; agents; 0]);
[t_sol, y_sol] = ode45(@subdep2_dif,0:dt:totalTime,[0 ; 0 ; 0]);

% R_ss = y_sol(end,1);
% Ep_ss = y_sol(end,2);               E_ss = y_sol(end,3);
% X_ss = y_sol(end,4);               

% disp(['Diff eq terminal R value   = ' num2str(y_sol(end,1))]);
% disp(['Diff eq terminal X value   = ' num2str(y_sol(end,3))]);

%% ** Graph deterministic results **
figure('Name','Substrate Depletion Time course','NumberTitle','off');
scatter(t_sol,y_sol(:,1),3,'.b');                               hold on;
scatter(t_sol,y_sol(:,2),3,'.g');
scatter(t_sol,y_sol(:,3),3,'.r');
% axis([0 t_sol(end) 0 agents]);          
xlabel('time');         ylabel('# Agents');                     
legend('R','Ep','X');                                           hold off;

% % ** Graph State Space **
% figure('Name','State Space','NumberTitle','off');               
% scatter(y_sol(:,3),y_sol(:,1),'+k');       % plot deterministic trajectory
% axis tight
% xlabel('X');                    ylabel('R');                  

%% Symbolic calculations
syms R Ep X positive;
dR_sym = + k_b1 * X + k_b2 * Ep * X - k_d * R;
dEp_sym = + k_r * (agents - Ep) * R / (Km_r + (agents - Ep)) - k_f * Ep / (Km_f + Ep);
dX_sym = + k_s * S - k_b1 * X - k_b2 * Ep * X;

% Numerical solution of steady-state values
SS = vpasolve([dR_sym,dEp_sym,dX_sym],[R,Ep,X]);  % Num solve for equns = 0
R_ss = double(SS.R);
Ep_ss = double(SS.Ep);
X_ss = double(SS.X);

% ** Nullclines **
% R_nc_sym = solve(dR_sym == 0,R);         % depends on Ep, X
% Ep_nc_sym = solve(dEp_sym == 0,Ep);      % 2 entries: +/-, depends on R
% X_nc_sym = solve(dX_sym == 0,X);         % depends on Ep, S
% 
% tempX = subs(X_nc_sym,Ep,Ep_nc_sym(1));
% tempR = subs(R_nc_sym,[Ep,X],[Ep_nc_sym(1),tempX]);
% R_ss = double(solve(tempR==R,R));
% Ep_ss = double(subs(Ep_nc_sym(1),R,R_ss));
% X_ss = double(subs(X_nc_sym,Ep,Ep_ss));


% ** Linearize and calculate eigenvalues **
fprintf(1,'Steady-States:\t\t\tEigenvalues\n');
Jac = jacobian([dR_sym,dEp_sym,dX_sym],[R,Ep,X]);

for w=1:size(R_ss,1)
    J = double(subs(Jac,[R,Ep,X],[R_ss(w),Ep_ss(w),X_ss(w)]));
    lambdas = eig(J);
    fprintf(['R=' num2str(R_ss(w)) ', Ep=' num2str(Ep_ss(w)) ', X=' num2str(X_ss(w)) ':\n']);
    disp([num2str(lambdas)]);
end
clear w;

%% Spectral Analysis
fs = 1/dt;                                          % sampling frequency
t_osc = totalTime/10:dt:totalTime;                          X_osc = y_sol(t_osc(1)/dt+1:end,3)';
Ep_osc = y_sol(t_osc(1)/dt+1:end,2)';               R_osc = y_sol(t_osc(1)/dt+1:end,1)';

% figure;     plot(t_osc,X_osc,t_osc,Yp_osc,t_osc,Rp_osc);
% xlabel = ('time (sec)');        ylabel('#');

Osc = [R_osc; Ep_osc; X_osc];
f_th = abs(imag(lambdas(2))) / (2*pi);              % Theoretical oscillation frequency (Hz)
period_th = 1 / f_th;                               % Theoretical period (sec)

figure('Name','Spectral Analysis: DE','NumberTitle','off','Position',[1 1 350 800]);          

% do FFT on R, Ep, X data
n = 2^(nextpow2(length(X_osc))+4);                      
f = fs / n * (0:n/2);       

species = {'\underbar{Species R}', '\underbar{Species Ep}', '\underbar{Species X}'}; 

Z = zeros(size(species,2),n);               Zn = zeros(size(species,2),n);

for b=1:size(species,2)
    Z(b,:) = fft(Osc(b,:),n);
    Zn(b,:) = abs(Z(b,:)/n);
    subplot(3,1,b);                                             hold on;
    title(species{b},'FontName','Times New Roman','FontSize',12,'Interpreter','LaTex');  
    plot(1000*f,Zn(b,1:n/2+1));                                   
    
%     axis([0.5*f_th 2*f_th 0 12*b]);
%     xlim([500*f_th 1500*f_th]);
    xlim([0 4]);                    % Adjust this when Im(lambdas) = 0 
    ylim([0 10*Zn(b,find(abs(f-f_th)<=0.0002 & abs(f-f_th)>=0.00005,1))]);


    tempy = [0 50 100];                     tempf = 1000*[f_th f_th f_th];           
    plot(tempf,tempy,'r');                  
    set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
    xlabel('f (mHz)');                      ylabel('Power');
end
  
hold off;                               clear b tempy tempf;
