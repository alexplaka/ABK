function[min] = Gettime()
%This function calculates the current time as the number of minutes that
%have elapsed since midnight.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

% format long;
currenttime = clock;
min = currenttime(4)*60 + currenttime(5) + currenttime(6)/60;