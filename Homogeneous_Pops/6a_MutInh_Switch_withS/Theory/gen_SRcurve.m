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

agents = 110;

S_array = 0:0.5:20;

A_ss = zeros(size(S_array,2),3);
B_ss = zeros(size(S_array,2),3);

for n=1:size(S_array,2)
    
    S = S_array(n);
    fprintf(1,['\nS = ' num2str(S) ':\n']);
    
    % *** Set up differential equations of system using symbolic variables ***
    syms A B k_sa_sym k_sb_sym positive;

    k_sa_sym = k_a * (K_b^n_b + alpha_b * B^n_b) / (K_b^n_b + B^n_b);
    k_sb_sym = k_b * (K_a^n_a + alpha_a * (A * 1/(1+(S/K_s)^n_s))^n_a) / ...
        (K_a^n_a + (A * 1/(1+(S/K_s)^n_s))^n_a);

    dA_sym = + k_sa_sym - k_da * A;
    dB_sym = + k_sb_sym - k_db * B;

    % Nullclines
    A_nc_sym = solve(dA_sym == 0,A,'MaxDegree',4,'Real',true);
    B_nc_sym = solve(dB_sym == 0,B,'MaxDegree',4,'Real',true);

    Jac = jacobian([dA_sym,dB_sym],[A,B]);                  % Calculate Jacobian matrix

    % Calculate nullclines, direction field
    A_nc = double(subs(A_nc_sym,B,0:agents));
    B_nc = double(subs(B_nc_sym,A,0:agents));

    if size(B_nc,1) > 1                % Nullcline values must be positive!
        for j=size(B_nc,1):-1:1
            if isempty(find(sign(B_nc(j,:)) == -1, 1,'first')) == false
                B_nc(j,:) = [];
            end
        end
    end

    [A_ss_int, B_ss_int] = intersections(A_nc,0:agents,0:agents,B_nc);

    % Clean up duplicate ss data
    for z=size(B_ss_int,1):-1:2
        temp = B_ss_int(z) - B_ss_int(z-1);
        if abs(temp) < 0.0001
            B_ss_int(z) = [];
            A_ss_int(z) = [];
        end
    end

%     fprintf(1,'\nFixed Points:\t\tEigenvalues\n');
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

    if isempty(A_ss_int) == 0
        A_ss(n,:) = A_ss_int(:);
        B_ss(n,:) = B_ss_int(:);
    else
        A_ss(n,:) = NaN;
        B_ss(n,:) = NaN;
    end
    
end

%% Remove duplicate FPs from monostable systems
for i=1:size(A_ss,1)
   if A_ss(i,1) == A_ss(i,2)
       A_ss(i,2:3) = NaN;
       B_ss(i,2:3) = NaN;
   end
end

%% Plot deterministic SR curve
fig1 = figure('Name','SR curve','NumberTitle','off');
set(fig1,'Position',[501 1 500 406]);                       hold on;

for h=1:size(B_ss,2)
    if mod(h,2) == 0
        p(h) = plot(S_array,B_ss(:,h),'--r','LineWidth',2,'DisplayName','Unstable');
    else
    p(h) = plot(S_array,B_ss(:,h),'-b','LineWidth',2,'DisplayName','Stable');
    end
end

axis([S_array(1) S_array(end) 0 max(max(B_ss))+10]);        hold off;
set(gca,'XMinorTick','on','YMinorTick','on','Box','off');
xlabel('N_S');         ylabel('N_R^*');         

leg1 = legend([p(1) p(2)]);
set(leg1,'Location','SouthEast');
set(leg1,'FontName','Times New Roman','FontSize',9,...
    'EdgeColor',[0.95 0.95 0.95]);

%% Save variables
save('SR_curve.mat');