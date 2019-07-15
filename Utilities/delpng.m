function [] = delpng()
% This function deletes PNG files from the current directory.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

list_fig = ls('*.png');             % Apply listing rule for PNG images
number = size(list_fig,1);

for n = 1:number
    filename = strcat(list_fig(n,:));
    delete(filename);
end
