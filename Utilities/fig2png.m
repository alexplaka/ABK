function [] = fig2png(FIGfile1)
% If called with one argument 'FIGfile1' (must be entered as a string), the
% file is converted to a PNG image. Alternatively, if no input is given,
% all FIG plots are converted to PNG images in a given directory (whatever 
% 'pwd' is when function is called).

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

if nargin == 1                          % Is a single FIG file specified?
    hgload(FIGfile1);
    print('-dpng',[FIGfile1 '.png']);
    close;
elseif nargin == 0                      % If no input is given, convert all
                                        % FIG files in pwd
    list_fig = ls('*.fig');             % Select FIG files
    number = size(list_fig,1);

    for n = 1:number
        filename = strcat(list_fig(n,:));
        hgload(filename);
        print('-dpng',[filename '.png']);
    end
    close all;
end
