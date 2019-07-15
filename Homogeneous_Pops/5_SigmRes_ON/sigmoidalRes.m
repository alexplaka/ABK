% Sigmoidal Response system (0th order ultrasensitivity) 
% (Michaelis-Menten kinetics for forward and reverse rxs)
% R <--> Rp             Rates: k_f, k_r
%     ^                 Michaelis-Menten constants: Km_f, Km_r
%     |
%     S                 S: Signal (reactant in forward reaction)
% Assume number of S molecules does NOT change.
% Considering individual agents.

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;             tic;
clc; 
%% Declare variables and functions
global k_f k_r Km_f Km_r S;


maxTime = 1000;                    % Maximum Simulation time (sec)
dt = 1/50;                         % Constant (fixed) time step increment (sec)
t_steps = maxTime / dt;

agents = 100;               disp(['Agents = ' num2str(agents)]);
k_f = 0.05;                             % basal forward rate constant
k_r = 1;                              % basal reverse rate constant
Km_f = 30;                              % MICROSCOPIC MM constant for forward rx
Km_r = 5;                              % MICROSCOPIC MM constant for reverse rx

S = 29;                            % Number of S molecules is constant
% *** Set up differential equations of system using symbolic variables ***
syms R Rp positive;
dR_sym =  - k_f * R * S / (Km_f + R)  +  k_r * Rp / (Km_r + Rp);
dRp_sym = + k_f * R * S / (Km_f + R)  -  k_r * Rp / (Km_r + Rp);

% % Nullclines
% R_nullcline_sym = solve(dR_sym == 0,R);
% Rp_nullcline_sym = solve(dRp_sym == 0,Rp);

%% Symbolic calculations
% Calculate theoretical steady state values for R, Rp
dR_sym_sub = subs(dR_sym,Rp,agents-R);       % Restriction: R+Rp=agents
R_ss_sym = solve(dR_sym_sub == 0,R,'MaxDegree',4,'Real',true);
R_ss = double(R_ss_sym);

% Only want positive steady state values and within right range (<agents)
if size(R_ss,1)>1
    R_ss_all = R_ss;
    temp = find(R_ss <= agents & R_ss >= 0,1);
    R_ss = R_ss_all(temp);
end
Rp_ss = agents - R_ss;
disp(['Theoretical  R_ss = ' num2str(R_ss)]);
disp(['Theoretical Rp_ss = ' num2str(Rp_ss)]);

% Linearize and calculate eigenvalues
Jac = jacobian([dR_sym,dRp_sym],[R,Rp]);
J = double(subs(Jac,[R,Rp],[R_ss,Rp_ss]));
TrJ = J(1,1) + J(2,2);
DetJ = det(J);
DiscrJ = TrJ^2 - 4*DetJ;
lambda_plus = (TrJ + sqrt(DiscrJ)) / 2;
lambda_minus = (TrJ - sqrt(DiscrJ)) / 2;        

%% Construct Rate curve
Rp_array = 0:agents;

dRpdt_degr = k_r .* Rp_array ./ (Km_r + Rp_array);
dRpdt_synt = k_f * S .* (agents - Rp_array) ./ (Km_f + (agents - Rp_array)); 

[Rp_ss_int, ssr] = intersections(Rp_array,dRpdt_degr,Rp_array,dRpdt_synt);
% Rp_ss_int: Steady-state value of Rp derived from intersections in Rate curve
% ssr: Steady-state rate of forward and back reactions
figure('Name','Rate Curve','NumberTitle','off');
plot(Rp_array,dRpdt_degr,'r',Rp_array,dRpdt_synt,'-k',Rp_ss_int,ssr,'ob');
xlabel('R_p');        ylabel('dR_p/dt');    
legend('degradation','synthesis','R_{p-ss}','Location','Best');
clear Rp_array;

%% Initialize - Preallocate memory for variables; ** Initial conditions **
t = 1;                              % Time counter variable
P_f = zeros(1,t_steps);             % For storing probability value at each time step
P_r = zeros(1,t_steps);             % For storing probability value at each time step
Tr = zeros(1,t_steps);              % Total sum of R agents in each time step	
Trp = zeros(1,t_steps);             % Total sum of Rp agents in each time step

% ******** Initial conditions - Number of R, Rp Agents ********
Tr(1) = floor(agents/2);       		Trp(1) = agents - Tr(1);					
% ************************************************************

tempR = zeros(1,agents);                          
% Put "1" where agents are "alive", then randomize the array
for c=1:Tr(1),                      tempR(c)=1;             end
tempR = RandArray(tempR);           % Randomize array
Rv = [tempR ; tempR];               % Initialize vector storing state of R agents     
Rpv = ~Rv;                          % Initialize vector storing state of Rp agents 
% Notes on Rv, Rpv:
% - Markov process, so only previous and current time steps needed --> 2 rows:          
% - Rv and Rpv are complementary (R+Rp=agents)
clear c tempR;

%% ABK Simulation

while t*dt <= maxTime               % Alt Condition: abs(Tr(t)-R_ss) > 0.5
    	
    P_f(t) = k_f * S / (Km_f + Tr(t)) * dt;  % P_dif forward rx
    P_r(t) = k_r / (Km_r + Trp(t)) * dt;     % P_dif reverse rx
    
    for i = 1:agents
        if Rv(1,i) == 1                      % if R agent is alive           
            if rand < P_f(t)                 % **check probability condition**
                Rv(2,i) = 0;                 % R agent "dies"
                Rpv(2,i) = 1;                % R is converted to Rp
            end
        elseif Rpv(1,i) == 1			     % if Rp agent is "alive"           
            if rand < P_r(t)                 % **check probability condition**
                Rpv(2,i) = 0;                % Rp agent "dies"
                Rv(2,i) = 1;                 % Rp is converted to R
            end
        end
    end
    Tr(t+1) = sum(Rv(2,:));                  Trp(t+1) = sum(Rpv(2,:));
    Rv(1,:) = Rv(2,:);                       Rpv(1,:) = Rpv(2,:); 
    t = t + 1;   
    
end

% % Remove unnecessary terminal 0's from arrays
% if t < t_steps
%     Tr = Tr(Tr~=0);                             w = size(Tr,2);                 
%     Trp = Trp(1:w);                             
%     P_f = P_f(1:w);               P_r = P_r(1:w);
% end

disp(['ABK sim terminal R value   = ' num2str(Tr(end))]);  
disp(['ABK sim terminal Rp value  = ' num2str(Trp(end))]);        
clear Rv Rpv;
finaltime = (t-1) * dt;

%% Solve differential equation (was waiting for variable 'finaltime')
[t_sol, y_sol] = ode45(@sigmoRes_dif,0:finaltime/100:finaltime,[Tr(1) ; Trp(1)]);
disp(['Diff eq terminal R value   = ' num2str(y_sol(end,1))]);
disp(['Diff eq terminal Rp value  = ' num2str(y_sol(end,2))]);

%% Graph ABK (stochastic) and deterministic results
figure('Name','Time course','NumberTitle','off');               hold on;
time = 0:dt:finaltime;
scatter(time,Tr,3,'oc');                                        
scatter(time,Trp,3,'om');
scatter(t_sol,y_sol(:,1),3,'.b');                  
scatter(t_sol,y_sol(:,2),3,'.r');                               hold off;
axis([0 t_sol(end) 0 agents]);          xlabel('time');         ylabel('# Agents');         
legend('R stoc','Rp stoc','R deter','Rp deter','Location','Best');

%% Finish
clear r w i j k temp;
toc
