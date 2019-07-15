function [s] = ElapsedTime(seconds);

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

et_sec = seconds;   et_min = et_sec / 60;  et_hr = et_min / 60;
if et_sec > 3600
%     fprintf(1,'\nElapsed time = %3.2f hrs \n',et_hr);
    s = sprintf('Elapsed time = %3.2f hrs \n',et_hr);
else
%     fprintf(1,'\nElapsed time = %3.2f min \n',et_min);
    s = sprintf('Elapsed time = %3.2f min \n',et_min);
end