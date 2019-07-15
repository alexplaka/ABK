%  A --> B
% Simulating 1st order kinetics. Considering individual agents.
% This script prepares the figure for inserting into the report. 

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

clear;       tic;
    
lambda = 0.8;                   k = log(2)/lambda;
totalTime = 8;             % Simulation time (sec)
dt = 1/10;                  % Constant (fixed) time step increment (sec)
agents = [500 100 25 10];   % size(agents,2) must be even!

if exist('dt') == 1         
    t_steps = totalTime / dt;
else
    t_steps = 1200;             % empirically chosen value (for lambda = 0.8)
end
figure('Position',[1 1 800 800]);

for a=1:size(agents,2)
    
    S = zeros(2,t_steps);               S(:,1) = agents(a);
    
    for b = 1:2                 % Do 2 runs of each simulation
    
        A = ones(t_steps,agents(a));          
        time = zeros(1,t_steps);
        P = zeros(1,t_steps);

        for t = 2:t_steps
    %         dt = 1 / (k*S(t-1));            % Variable time step increment
        %     dt = exprnd(1/(k*S(t-1)));   % Exponentially-distributed time step increment
            time(t) = time(t-1) + dt;      
            for i = 1:size(A,2)
                if A(t,i) == 1                      % if agent is still alive
                    P(t) = 1 - 2^(-dt/lambda);      % P from integrated rate law
        %             P(t) = k*dt;                  % P from differential rate law
                    if rand < P(t)                  % **check probability condition**
                        A(t,i) = 0;                 % agent dies
                        for w = t+1:t_steps
                            A(w,i) = 0;             % mark agent dead for rest of time
                        end
                    end
                end
            end
            S(b,t) = sum(A(t,:));
        end
        clear t w i;
    end
    
    subplot(size(agents,2)/2,size(agents,2)/2,a);
%     k_est = FitCurve_exp_FixA_g(time',S',agents(a));
    plot(time,S(1,:),'ob','MarkerSize',3);                  hold on;
    plot(time,S(2,:),'+c','MarkerSize',3);
    x = 0:0.01:totalTime;        y = agents(a) * exp(-k*x);
    plot(x,y,'r');
    xlabel('t (sec)');           ylabel('N_A'); 
    hleg = legend('Run 1','Run 2',['N_A(t) =' num2str(agents(a)) '\cdote^{-k\cdot\Deltat}']);
    set(hleg,'FontName','Times New Roman','FontSize',8','Location','NorthEast');
    axis('tight');               hold off;

end

toc
