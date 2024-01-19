%    k_da           k_db
%    <--  A _  _ B  -->                 % For Report: B = R
%         ^  \/  ^
%     k_a |  /\  | k_b
%         | / ~\ |
%         |~  | ~|
%             |
%             S
% A and B are synthesized (0th order constants k_a, k_b respectively)
% A and B are degraded (k_da, k_db)
% A and B inhibit each other's synthesis (MUTUAL inhibition, alpha < 1):
% B influences the rate of the synthesis of A
% A influences the rate of the synthesis of B
% S inhibits the inhibition of B by A (assume complete repression)
% Assume S does not change.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

%% Declare variables and functions
clear;                          clc; 

global k_da k_db agents;
global k_a alpha_a K_b n_b;
global k_b alpha_b K_a n_a;
global K_s n_s S;

maxTime = 2000;                        % Maximum Simulation time (sec)
agents = 150;

k_a = 1;                            % MICROSCOPIC basal rate of A synthesis (0th order)
k_b = 1;                            % MICROSCOPIC basal rate of B synthesis (0th order)
k_da = 0.01;                        % degradation rate constant (1st order) for A
k_db = 0.01;                        % degradation rate constant (1st order) for B
% Populations producing half-maximal regulatory effect
K_b = 50;                   K_a =  50;
K_s = 50;

% *** Feedback parameters ***
% Degree of activation: >1 activator, <1 repressor, =1 no regulation
alpha_a = 0;                        % Effect of A on synthesis of B; = 0 means Complete Repression
alpha_b = 0;                        % Effect of B on synthesis of A; = 0 means Complete Repression
% Hill Coefficient: measure of cooperativity
n_b = 3;                      n_a = 3;
n_s = 1;

% *** Initial condition ***
Ai = 10;                    Bi = 10;

S = 30;

% *** Set up differential equations of system using symbolic variables ***
syms A B k_sa_sym k_sb_sym positive;

k_sa_sym = k_a * (K_b^n_b + alpha_b * B^n_b) / (K_b^n_b + B^n_b);
k_sb_sym = k_b * (K_a^n_a + alpha_a * (A * 1/(1+(S/K_s)^n_s))^n_a) / (K_a^n_a + (A * 1/(1+(S/K_s)^n_s))^n_a);

dA_sym = + k_sa_sym - k_da * A;
dB_sym = + k_sb_sym - k_db * B;

% Nullclines
A_nc_sym = solve(dA_sym == 0,A,'MaxDegree',4,'Real',true);
B_nc_sym = solve(dB_sym == 0,B,'MaxDegree',4,'Real',true);

Jac = jacobian([dA_sym,dB_sym],[A,B]);                  % Calculate Jacobian matrix

%% Calculate nullclines, direction field
A_nc = double(subs(A_nc_sym,B,0:agents));
B_nc = double(subs(B_nc_sym,A,0:agents));

if size(B_nc,1) > 1                % Nullcline values must be positive!
    for j=size(B_nc,1):-1:1
        if isempty(find(sign(B_nc(j,:)) == -1, 1,'first')) == false
            B_nc(j,:) = [];
        end
    end
end

[Am,Bm] = meshgrid(1:agents/15:agents);     % Mesh-grid values for constructing state space
k_sa_num = subs(k_sa_sym,[A,B],{Am,Bm});
k_sb_num = subs(k_sb_sym,[A,B],{Am,Bm});
dA_num = double(subs(dA_sym,[k_sa_sym,k_sb_sym,A,B],{k_sa_num,k_sb_num,Am,Bm}));
dB_num = double(subs(dB_sym,[k_sa_sym,k_sb_sym,A,B],{k_sa_num,k_sb_num,Am,Bm}));
% r = ( dA_num.^2 + dB_num.^2 ).^0.5;
r = 1;
Ua = dA_num./r;          Ub = dB_num./r;
V = sqrt(Ua.^2 + Ub.^2);                % Magnitude of velocity vector

[A_ss_int, B_ss_int] = intersections(A_nc,0:agents,0:agents,B_nc);

% Clean up duplicate ss data
for z=size(B_ss_int,1):-1:2
    temp = B_ss_int(z) - B_ss_int(z-1);
    if abs(temp) < 0.0001
        B_ss_int(z) = [];
        A_ss_int(z) = [];
    end
end

fprintf(1,'\nFixed Points:\t\tEigenvalues\n');
for w=1:size(A_ss_int,1)
    J = double(subs(Jac,[A,B],[A_ss_int(w),B_ss_int(w)]));
    TrJ = J(1,1) + J(2,2);
    DetJ = det(J);
    DiscrJ = TrJ^2 - 4*DetJ;
    lambda_plus(w) = (TrJ + sqrt(DiscrJ)) / 2;
    lambda_minus(w) = (TrJ - sqrt(DiscrJ)) / 2;   
    fprintf(1,['A= %4.2f , B= %4.2f :\t%+6.4f , %+6.4f \n'],...
        A_ss_int(w),B_ss_int(w),lambda_plus(w),lambda_minus(w));
end

%% Solve differential equation (was waiting for variable 'finaltime')
[t_sol, y_sol] = ode45(@mutual_withS_dif,0:maxTime/100:maxTime,[Ai ; Bi]);

%% Plot State Space: nullclines, direction field, trajectories, separatrix, fixed points
fig1 = figure('Name','State Space','NumberTitle','off');                hold on;
set(fig1,'Position',[1 1 500 406]); 

p_ncA = plot(A_nc,0:agents,'','Color',[1 0.3 0.2],...
    'LineWidth',2,'DisplayName','N_A nullcline');
p_ncR = plot(0:agents,B_nc,'','Color',[0.4 0.4 1],...
    'LineWidth',2,'DisplayName','N_R nullcline');      % B = R

streakarrow(Am,Bm,Ua,Ub,0.7,1);             % Direction field
% colorbar vert;             h=colorbar;    set(h,'ylim',[min(Vmag(:)) max(Vmag(:))]);
% quiver(Am,Bm,U,V,'g');                    % Direction field

% Plot separatrix
% p_s = plot(0:agents,0:agents,'--','Color',[1 0.75 0],...
%     'LineWidth',1.5,'DisplayName','Separatrix');

% Plot DE trajectory
% p_de  = plot(y_sol(:,1),y_sol(:,2),'+k','DisplayName','DE  Trajectory','MarkerSize',8);  

if size(A_ss_int,1)==3
    for z=1:3
        if mod(z,2)==0              % Unstable fixed point
            p(z) = plot(A_ss_int(z),B_ss_int(z),'oc','MarkerSize',8,...
                'MarkerFaceColor','r','DisplayName','Unstable FP');
        else                        % Two Stable fixed points
            p(z) = plot(A_ss_int(z),B_ss_int(z),'oc','MarkerSize',8,...
                'MarkerFaceColor','b','DisplayName','Stable FP');
        end
    end
else
    p = plot(A_ss_int,B_ss_int,'oc','MarkerSize',8,...
        'MarkerFaceColor','b','DisplayName','Stable FP');          % Sole stable fixed point
end

xlabel('N_A');                    ylabel('N_R');  
title(['N_S = ' num2str(S)],'FontName','Times New Roman','FontSize',11);
axis([0 110 0 110]);                                             hold off;  
set(gca,'XMinorTick','on','YMinorTick','on','Box','off','DataAspectRatio',[1 1 1]);

leg1 = legend([p_ncA p_ncR p(:)]);    
set(leg1,'Location','NorthEast');
% set(leg1,'Position',[0.7 0.71 0.29 0.262]);
set(leg1,'FontName','Times New Roman','FontSize',8,...
    'EdgeColor',[0.95 0.95 0.95]);
