function [] = rfig2png()
% Recursive fig2png. Subdirectories of current folder are scanned for FIG
% files which are converted to PNG images.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

fig2png;                            % Do conversion: fig2png in pwd
list = ls;
number = size(list,1);

for n = 3:number                    % Starts from 3 such that '.' and
                                    % '..' in 'ls' output are ignored
    filename = strcat(list(n,:));
    if isdir(filename) == true      % Is it a directory?
        cd(filename);               % If yes, change to it
        rfig2png;                   % Recursive call to rfig2png
        cd ..;                      % Return to parent directory
    end
end