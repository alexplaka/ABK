%    R  <--->  Rp            Rates: k_f, k_r (Forward rx activated by X)
%          ^    |            Michaelis-Menten constants: Km_f, Km_r
%         /     |
%         |     \
%    -->  X    -->                      Rp: Response
%   k_b   ^    k_d1/k_d2
%         | k_s
%         S                             S: Signal

% Simulating negative feedback 2-component process. 
% k_b: X synthesis, 0th order process
% k_s: X synthesis, 1st order process wrt S (but S is not consummed)
% k_s: ALTERNATIVE IMPLEMENTATION: X synthesis, 0th order process, UPregulated by S
% k_d1: X degradation, 1st order process wrt X
% k_d2: X degradation, 2nd order process wrt X, Rp (but Rp is not consummed)
% k_f: Rp synthesis, MM process wrt X (but X is not consummed)
% k_r: R synthesis, MM process 
% Km_f, Km_r: MM constants for forward and reverse rxs

% Note: Assume number of S molecules/agents is NOT changing (eg, S is DNA molecules)

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;  % clc;

global agents k_s k_b k_d1 k_d2 k_f k_r Km_f Km_r S;

agents = 100;
k_b = 0.0;                        % 0th order X synthesis rate
k_d1 = 0.00;                       % 1st order degradation rate (units: 1/sec)
k_d2 = 0.005;                       % 2nd order degradation rate wrt X, Rp (units: 1/sec)
k_s = 0.10;                         % 0th order R synthesis rate UPregulated by S
k_f = 0.005;                           % basal forward rate (1st order)
k_r = 0.5;                        % basal reverse rate (1st order)
Km_f = 30;                          % MICROSCOPIC Michaelis-Menten constant for forward rx
Km_r = 1;                           % MICROSCOPIC Michaelis-Menten constant for reverse rx
S = 30;                            % Assume number of S molecules/agents is NOT changing

% Solve DE
[t_sol, y_sol] = ode23(@negFb_dif,0:2500,[0 ; 0]);   
% t=0:2000 sec, to ensure that steady-states have been reached.

X_ss = y_sol(end,1);                Rp_ss = y_sol(end,2);

% disp(['Diff eq terminal X value   = ' num2str(y_sol(end,1))]);
% disp(['Diff eq terminal Rp value   = ' num2str(y_sol(end,2))]);

plot_tmax = 2000;
% Graph deterministic results
figure('Name','Negative Feedback Time course','NumberTitle','off');     hold on;
scatter(t_sol(1:plot_tmax),y_sol(1:plot_tmax,1),3,'.b');                               
scatter(t_sol(1:plot_tmax),y_sol(1:plot_tmax,2),3,'.r');
% axis([0 t_sol(end) 0 agents]);          

xlabel('t (sec)');         ylabel('N');                                 hold off;

% ** Symbolic calculations - Linearize and calculate eigenvalues **
syms X Rp positive;
dX_sym = + k_b + k_s * S - k_d1 * X - k_d2 * X * Rp;
dRp_sym = + k_f * (agents - Rp) * X / (Km_f + (agents-Rp)) - k_r * Rp / (Km_r + Rp);

fprintf(1,'Steady-States:\t\t\tEigenvalues\n');
Jac = jacobian([dX_sym,dRp_sym],[X,Rp]);                  % Calculate Jacobian matrix

for w=1:size(X_ss,1)
    J = double(subs(Jac,[X,Rp],[X_ss(w),Rp_ss(w)]));
    TrJ = J(1,1) + J(2,2);
    DetJ = det(J);
    DiscrJ = TrJ^2 - 4*DetJ;
    lambda_plus(w) = (TrJ + sqrt(DiscrJ)) / 2;
    lambda_minus(w) = (TrJ - sqrt(DiscrJ)) / 2;   
    fprintf(1,['X=' num2str(X_ss(w)) ', Rp=' num2str(Rp_ss(w)) ':\t' ...
        num2str(lambda_plus(w)) ', ' num2str(lambda_minus(w)) '\n']);
end
clear w;
