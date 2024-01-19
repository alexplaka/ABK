% Hyperbolic Response system:
% R <--> Rp             Rates: k_f, k_r
%     ^
%     |
%     S                 S: Signal (promotes forward reaction)
% Assume number of S molecules does NOT change.
% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear;             tic;
clc; 
rng(0);
%% Declare variables and functions
global k_f k_r S;

maxTime = 50;                          % Maximum Simulation time (sec)
dt = 1/50;                              % Constant (fixed) time step increment (sec)
t_steps = maxTime / dt;

agents = 100;          	disp(['Agents = ' num2str(agents)]);
k_f = 0.01;                            % basal forward rate (1st order)
k_r = 0.1;                             % basal reverse rate (1st order)
S = 1;                                 % Number of S molecules is constant

%% *** Set up differential equations of system using symbolic variables ***
% *** (requires Symbolic Math Toolbox) ***

% syms R Rp positive;
% dR_sym =  -k_f * R * S + k_r * Rp;
% dRp_sym =  k_f * R * S - k_r * Rp;
% 
% % Nullclines
% R_nullcline_sym = solve(dR_sym == 0,R);
% Rp_nullcline_sym = solve(dRp_sym == 0,Rp);
% 
% % ** Symbolic calculations **
% % Calculate theoretical steady state values for R, Rp
% dR_sym_sub = subs(dR_sym,Rp,agents-R);       % Restriction: R+Rp=agents
% R_ss_sym = solve(dR_sym_sub == 0,R,'MaxDegree',4,'Real',true);
% R_ss = double(R_ss_sym);
% % Pen & paper sol'n:  R_ss = agents / ( (k_f/k_r)*S + 1 );       Rp_ss = agents - R_ss;
% 
% % Only want positive steady state values
% if size(R_ss,1)>1
%     R_ss_all = R_ss;
%     temp = find(R_ss <= agents,1);
%     R_ss = R_ss_all(temp);
% end
% Rp_ss = agents - R_ss;
% disp(['Theoretical  R_ss = ' num2str(R_ss)]);
% disp(['Theoretical Rp_ss = ' num2str(Rp_ss)]);
% 
% % Linearize and calculate eigenvalues
% Jac = jacobian([dR_sym,dRp_sym],[R,Rp]);
% J = double(subs(Jac,[R,Rp],[R_ss,Rp_ss]));
% TrJ = J(1,1) + J(2,2);
% DetJ = det(J);
% DiscrJ = TrJ^2 - 4*DetJ;
% lambda_plus = (TrJ + sqrt(DiscrJ)) / 2;
% lambda_minus = (TrJ - sqrt(DiscrJ)) / 2;        

%% Simple formula for Rp_ss
Rp_ss = agents * S / ( k_r/k_f + S );

%% Construct Rate curve
Rp_array = 0:agents;

dRpdt_degr = k_r .* Rp_array;
dRpdt_synt = k_f * S .* (agents - Rp_array);  

[Rp_ss_int, ssr] = intersections(Rp_array,dRpdt_degr,Rp_array,dRpdt_synt,0);
% Rp_int: Steady-state value of Rp derived from intersections in Rate curve
% ssr: Steady-state rate of forward and back reactions
figure('Name','Rate Curve','NumberTitle','off');
plot(Rp_array,dRpdt_degr,'r',Rp_array,dRpdt_synt,'-k',Rp_ss_int,ssr,'ob');
xlabel('R_p');        ylabel('dR_p/dt');    
legend('degradation','synthesis','R_{p-ss}','Location','NorthEast');
clear Rp_array;

%% Initialize - Preallocate memory for variables; ** Initial conditions **
t = 1;                              % Time counter variable
Tr = zeros(1,t_steps);              % Total sum of R agents in each time step	
Trp = zeros(1,t_steps);             % Total sum of Rp agents in each time step

% ******** Initial conditions - Number of R, Rp Agents ********
Tr(1) = floor(agents/4);      		Trp(1) = agents - Tr(1);					
% ************************************************************

tempR = zeros(1,agents);                          
% Put "1" where agenst are "alive", then randomize the array
for c=1:Tr(1),                      tempR(c)=1;             end
tempR = RandArray(tempR);           % Randomize array
Rv = [tempR ; tempR];               % Initialize vector storing state of R agents     
Rpv = ~Rv;                          % Initialize vector storing state of Rp agents 
% Notes on Rv, Rpv:
% - Markov process, so only previous and current time steps needed --> 2 rows:          
% - Rv and Rpv are complementary (R+Rp=agents)
clear c tempR;

%% ABK Simulation

P_f = 1 - exp(-k_f * S * dt);    % P_ber for forward reaction: 1st order            
%     P_f = k_f * S * dt;              % P_dif for forward reaction: 1st order
P_r = 1 - exp(-k_r * dt);        % P_ber for reverse reaction: 1st order            
%     P_r = k_r * dt;                  % P_dif for reverse reactio: 1st order

% while abs(Trp(t)-Rp_ss) > 0.5 && t*dt <= maxTime
while t*dt <= maxTime
    
%     dt = 1 / (ko*Sa(t-1));            % Variable time step increment ; OLD
%     dt = exprnd(1/(k*S(t-1)));   % Exponentially-distributed time step increment ; OLD
    
    for i = 1:agents
        if Rv(1,i) == 1                     % if R agent is alive
            if rand < P_f                  % **check probability condition**
                Rv(2,i) = 0;                 % R agent "dies"
                Rpv(2,i) = 1;                 % R is converted to Rp
            end
        elseif Rpv(1,i) == 1					% if Rp agent is "alive"
            if rand < P_r                  % **check probability condition**
                Rpv(2,i) = 0;                 % Rp agent "dies"
                Rv(2,i) = 1;                 % Rp is converted to R
            end
        end
    end
    Tr(t+1) = sum(Rv(2,:));                  Trp(t+1) = sum(Rpv(2,:));
    Rv(1,:) = Rv(2,:);                       Rpv(1,:) = Rpv(2,:); 
    t = t + 1;   
end

% Remove unnecessary terminal 0's from arrays
Tr = Tr(Tr~=0);                             w = size(Tr,2);                 
Trp = Trp(1:w);                                      
disp(['ABK sim terminal R value   = ' num2str(Tr(end))]);  
disp(['ABK sim terminal Rp value  = ' num2str(Trp(end))]);        
clear Rv Rpv;
finaltime = (t-1) * dt;

%% Solve differential equation (was waiting for variable 'finaltime')
[t_sol, y_sol] = ode45(@hyperRes_dif,0:finaltime/100:finaltime+(finaltime/20),[Tr(1) ; Trp(1)]);
disp(['Diff eq terminal R  value   = ' num2str(y_sol(end,1))]);
disp(['Diff eq terminal Rp value   = ' num2str(y_sol(end,2))]);
% Note: extended time domain of solution by finaltime/20 to get closer to steady-state values

%% Graph ABK (stochastic) and deterministic results
figure('Name','Time course','NumberTitle','off');
time = 0:dt:finaltime;
scatter(time,Tr,3,'oc');                                        hold on;
scatter(time,Trp,3,'om');
scatter(t_sol,y_sol(:,1),3,'.b');                  
scatter(t_sol,y_sol(:,2),3,'.r');                               hold off;
axis([0 t_sol(end) 0 agents]);          xlabel('time');         ylabel('# Agents');         
legend('R stoc','Rp stoc','R deter','Rp deter');

%% Plot State Space: nullclines, direction field, trajectories

% R_nullcline  = double(subs(R_nullcline_sym,Rp,0:agents));
% Rp_nullcline = double(subs(Rp_nullcline_sym,R,0:agents));
% 
% if size(R_nullcline,1)>1                % Nullcline values must be positive!
%     for j=size(R_nullcline,1):-1:1
%         if isempty(find(sign(R_nullcline(j,:)) == -1, 1,'first')) == false
%             R_nullcline(j,:) = [];
%         end
%     end
% end
% 
% if size(Rp_nullcline,1)>1               % Nullcline values must be positive!
%     for k=size(Rp_nullcline,1):-1:1
%         if isempty(find(sign(Rp_nullcline(k,:)) == -1, 1,'first')) == false
%             Rp_nullcline(k,:) = [];
%         end
%     end
% end
% 
% [Rm,Rpm] = meshgrid(0:agents/20:agents);     % Mesh-grid values for constructing state space
% dR_num  = double(subs(dR_sym,[R,Rp],{Rm,Rpm}));
% dRp_num = double(subs(dRp_sym,[R,Rp],{Rm,Rpm}));
% % r = ( dR_num.^2 + dRp_num.^2 ).^0.5;
% r=1;
% 
% figure('Name','State Space','NumberTitle','off');               hold on;
% plot(R_nullcline,0:agents,'b-.',0:agents,Rp_nullcline,'m-.','LineWidth',1);
% quiver(Rm,Rpm,dR_num./r,dRp_num./r,'g');                          axis([0 agents 0 agents]);                  
% scatter(y_sol(:,1),y_sol(:,2),'+k');                    % plot deterministic trajectory
% scatter(Tr,Trp,3,'.r');                                  % plot stochastic trajectory
% xlabel('R');                    ylabel('R_p');                  hold off; 
% legend('R nullcline','R_p nullcline','Direction field','Deterministic','ABK stochastic');

%% Finish
clear r w i j k temp;
toc
