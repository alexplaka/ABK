%% Calculate and plot SR curve.

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear; 

agents = 100;          	  % disp(['Agents = ' num2str(agents)]);

k_f = 0.05;                             % basal forward rate constant
Km_f = 30;                              % MICROSCOPIC MM constant for forward rx
k_r = 1;                              % basal reverse rate constant
Km_r = 5;                               % MICROSCOPIC MM constant for reverse rx

S_array = 1:60;

% Symbolic calcs; Plot SR curve
syms R Rp S positive;
dR_sym =  - k_f * R * S / (Km_f + R)  +  k_r * Rp / (Km_r + Rp);
temp = subs(dR_sym,R,agents-Rp);
Rp_sym = solve(temp==0,Rp);                 % result has 2 entries: -/+

Rp_ss = zeros(1,size(S_array,2));

for n=1:size(S_array,2)
    if n==k_r/k_f,     Rp_ss(n)=NaN;       continue;       end
    Rp_ss(n) = double(subs(Rp_sym(2),S,S_array(n)));
end

% ** Renaming parameters for simplicity **
T = agents;
a = k_f;
b = Km_f;
c = k_r;
d = Km_r;
% ****************************************

% Point of steepest ascent at (S_50, T/2)
S_50 = c/a * (T+2*b) / (T+2*d)
%% Plot
fig1 = figure;
set(fig1,'Position',[1 1 500 406]);                     hold on;       
plot(Rp_ss);            
plot(S_50,T/2,'or');

% - Indeed, plotting Rp_ss vs S produces a sigmoidal curve, as previously predicted. 
% - Note that this simple code can predict the theoretical dependence of a species' 
% steady-state value as a function of a parameter of choice (in this case, S).

%% Estimate greatest slope of ascent (at Rp = T/2) and compare to theoretical expectation
upper_end = find(abs(Rp_ss - (T/2 + T/20)) < 5);
lower_end = find(abs(Rp_ss - (T/2 - T/20)) < 5);

if size(upper_end,2) > 1
    upper_end = upper_end(1);
end

if size(lower_end,2) > 1
    lower_end = lower_end(end);
end

gslope_est = (Rp_ss(upper_end) - Rp_ss(lower_end)) / (S_array(upper_end) - S_array(lower_end))

gslope = a/4/c *T*(T+2*d)^2 / ( b*(T+2*d) + d*(T+2*b) )
