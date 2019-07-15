function [] = rdelpng()
% Recursive delpng. Subdirectories of current folder are scanned for PNG
% images which are deleted.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

delpng;                             % Do conversion: fig2png in pwd
list = ls; 
number = size(list,1);

for n = 3:number                    % Starts from 3 such that '.' and
                                    % '..' in 'ls' output are ignored
    filename = strcat(list(n,:));
    if isdir(filename) == true      % Is it a directory?
        cd(filename);               % If yes, change to it
        rdelpng;                    % Recursive call to rdelpng
        cd ..;                      % Return to parent directory
    end
end