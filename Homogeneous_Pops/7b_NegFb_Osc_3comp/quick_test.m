%          R <---->  Rp             Rates: k_f, k_r (Forward rx activated by Yp)
%               ^     |             Michaelis-Menten constants: Km_f, Km_r
%              /      |
%             |       |
%    Y <---> Yp       |             Rates: k_f1, k_r1 (Forward rx activated by X)
%         ^           |             Michaelis-Menten constants: Km_f1, Km_r1
%        /            |
%        |            \
%   -->  X  ---------------->                      Rp: Response
%   k_b  ^       k_d1/k_d2
%        | 
%        | k_s
%        S                             S: Signal

% Simulating negative feedback process (3-component system: X, Y/Yp, R/Rp) 
% k_b: X synthesis, 0th order process
% k_s: X synthesis, 1st order process wrt S (but S is not consummed)
% k_s: ALTERNATIVE IMPLEMENTATION: X synthesis, 0th order process, UPregulated by S
% k_d1: X degradation, 1st order process wrt X
% k_d2: X degradation, 2nd order process wrt X, Rp (but Rp is not consummed)
% k_f: Rp synthesis, MM process wrt Yp (but Yp is not consummed)
% k_r: R synthesis, MM process 
% Km_f, Km_r: MM constants for forward and reverse rxs R <==> Rp
% k_f1: Yp synthesis, MM process wrt X (but X is not consummed)
% k_r1: Y synthesis, MM process 
% Km_f1, Km_r1: MM constants for forward and reverse rxs Y <==> Yp

% Note: Assume number of S molecules/agents is NOT changing (eg, S is DNA molecules)

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear; clc;
global agents k_s k_b k_d1 k_d2 k_f k_r Km_f Km_r k_f1 k_r1 Km_f1 Km_r1 S;

agents = 250;
k_b = 0;                            % 0th order X synthesis rate
k_d1 = 0.00;                       % 1st order degradation rate (units: 1/sec)
k_d2 = 0.001;                       % 2nd order degradation rate wrt X, Rp (units: 1/sec)
k_s = 0.40;                         % 0th order R synthesis rate UPregulated by S

k_f = 0.01;                         % basal forward rate R-->Rp
k_r = 1;                            % basal reverse rate R<--Rp
Km_f = 10;                          % MICROSCOPIC Michaelis-Menten constant for forward rx R-->Rp
Km_r = 10;                          % MICROSCOPIC Michaelis-Menten constant for reverse rx R<--Rp

k_f1 = 0.01;                        % basal forward rate Y-->Yp
k_r1 = 1;                           % basal reverse rate Y<--Yp
Km_f1 = 10;                         % MICROSCOPIC Michaelis-Menten constant for forward rx Y-->Yp
Km_r1 = 10;                         % MICROSCOPIC Michaelis-Menten constant for reverse rx Y<--Yp

S = 30;                             % Assume number of S molecules/agents is NOT changing

% Solve DE: Choose between the 2 implementations below (3 vs 5 species)
% [t_sol, y_sol] = ode15s(@negFb3_dif,0:2500,[0 ; 100; 100]);
[t_sol, y_sol] = ode45(@negFb3_dif2,0:10000,[0 ; 150; 100; 150; 100]);
% *** Note: Only Implementation #2 works! ***

% disp(['Diff eq terminal X value   = ' num2str(y_sol(end,1))]);
% disp(['Diff eq terminal Rp value   = ' num2str(y_sol(end,3))]);

%% Graph deterministic results
figure('Name','Negative Feedback Time course','NumberTitle','off');
scatter(t_sol,y_sol(:,1),3,'.b');                               hold on;
scatter(t_sol,y_sol(:,3),3,'.g');
scatter(t_sol,y_sol(:,5),3,'.r');
% axis([0 t_sol(end) 0 agents]);   
axis tight;
xlabel('time');         ylabel('# Agents');                     
legend('X','Yp','Rp');                                          hold off;

%% Symbolic calculations - Linearize and calculate eigenvalues
syms X Y Yp R Rp positive;
dX_sym  = + k_b + k_s * S - k_d1 * X - k_d2 * X * Rp;
dY_sym  = - k_f1 * Y * X / (Km_f1 + Y) + k_r1 * Yp / (Km_r1 + Yp);
dYp_sym = + k_f1 * Y * X / (Km_f1 + Y) - k_r1 * Yp / (Km_r1 + Yp);
dR_sym  = - k_f * R * Yp / (Km_f + R) + k_r * Rp / (Km_r + Rp);
dRp_sym = + k_f * R * Yp / (Km_f + R) - k_r * Rp / (Km_r + Rp);
% Nullclines
% X_nc_sym = solve(dX_sym == 0,X,'MaxDegree',4,'Real',true);
% Yp_nc_sym = solve(dYp_sym == 0,Yp,'MaxDegree',4,'Real',true);
% Rp_nc_sym = solve(dRp_sym == 0,Rp,'MaxDegree',4,'Real',true);

%% 
X_ss = y_sol(end,1);
Y_ss = y_sol(end,2);               Yp_ss = y_sol(end,3);
R_ss = y_sol(end,4);               Rp_ss = y_sol(end,5);

%%
fprintf(1,'Steady-States:\t\t\tEigenvalues\n');
Jac = jacobian([dX_sym,dY_sym,dYp_sym,dR_sym,dRp_sym],[X,Y,Yp,R,Rp]);

for w=1:size(X_ss,1)
    J = double(subs(Jac,[X,Y,Yp,R,Rp],[X_ss(w),Y_ss(w),Yp_ss(w),R_ss(w),Rp_ss(w)]));
    lambdas = eig(J);
    fprintf(['X=' num2str(X_ss(w)) ', Yp=' num2str(Yp_ss(w)) ', Rp=' num2str(Rp_ss(w)) ':\n']);
    disp([num2str(lambdas)]);
end
clear w;
