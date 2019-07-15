function[c_time] = CreationTime()
% This function creates a string in the following format:
% ******************************************************
% month-day-'minutes since midnight (day start=0min)'_
% ******************************************************

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

currenttime = clock;
month = currenttime(2);
day = currenttime(3);
min = round(Gettime);
c_time = [num2str(month) '-' int2str(day) '-' int2str(min) '_'];