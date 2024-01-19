% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

clear; 

agents = 100;          	  % disp(['Agents = ' num2str(agents)]);

k_f = 1;                             % basal forward rate constant
Km_f = 5;                              % MICROSCOPIC MM constant for forward rx
k_r = 0.05;                              % basal reverse rate constant
Km_r = 30;                               % MICROSCOPIC MM constant for reverse rx

S_array = 1:60;

% Symbolic calcs; Plot SR curve
syms R Rp S positive;
dRp_sym =  + k_f * R / (Km_f + R)  -  k_r * S * Rp / (Km_r + Rp);
temp = subs(dRp_sym,R,agents-Rp);
Rp_sym = solve(temp==0,Rp,'MaxDegree',4,'Real',true);                 % result has 2 entries 

Rp_ss = zeros(1,size(S_array,2));

for n=1:size(S_array,2)
    if S_array(n)==k_f/k_r,     Rp_ss(n)=NaN;       continue;       end
    Rp_ss(n) = double(subs(Rp_sym(1),S,S_array(n)));
end

% ** Renaming parameters for simplicity **
T = agents;
a = k_f;
b = Km_f;
c = k_r;
d = Km_r;
% ****************************************

% Point of steepest ascent at (S_50, T/2)
S_50 = ( c/a * (T+2*b) / (T+2*d) )^-1
% Plot
fig1 = figure;
set(fig1,'Position',[1 1 500 406]);                     hold on;       
plot(Rp_ss);            
plot(S_50,T/2,'or');

% - Indeed, plotting Rp_ss vs S produces a sigmoidal curve, as previously predicted. 

gslope = c/4/a .* T.* (T+2*b).^2 ./ ( b .* (T+2*d) + d .* (T+2*b) )