% This script calculates pi by comparing the random 'hits' inside the unit
% circle centered at the origin and the inscribing square. The ratio of a
% circle's area to that of a square is pi/4.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

clear; clc;
%% Initialize
trials = 10000;
point = zeros(trials,2);
count = 0;      % Number of 'hits' inside unit circle

%% Main
for h = 1:trials
    
    sign = rand(1,2);
    if sign(1) < 0.5    sign_x = +1;
    else sign_x = -1;   end;
    if sign(2) < 0.5    sign_y = +1;
    else sign_y = -1;   end;

    point(h,:) = [sign_x sign_y] .* rand(1,2);

    if point(h,1)^2+point(h,2)^2 <= 1
        count = count + 1;
    end
    
end

pi_calc = count * 4 / trials
perc_error = 100 * abs(pi_calc - pi)/pi

%% Plot
figure; hold on;
scatter(point(:,1),point(:,2),1,'r');
ezplot('x^2+y^2=1');

%% Finish
clear h sign sign_x sign_y;