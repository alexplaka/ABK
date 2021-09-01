function [] = RateCurve(agents)

global k_b k_d k_s S;

R = 0:agents;

dRdt_degr = k_d .* R;
dRdt_synt = k_b + k_s * S * ones(1,size(dRdt_degr,2));  

[R_ss, r] = intersections(R,dRdt_degr,R,dRdt_synt);
figure('Name','Rate Curve','NumberTitle','off');
plot(R,dRdt_degr,'r',R,dRdt_synt,'-k',R_ss,r,'ob');
xlabel('R');        ylabel('dR/dt');    
legend('degradation','synthesis','R_{ss}','Location','NorthWest');