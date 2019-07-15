function rA = RandMatrix(A)
% This function randomizes the order of an matrix's contents.
% The Matlab 'sort' function only sorts along rows and columns separately,
% so I have adopted a sequential sorting strategy, where a random matrix is
% first sorted along its columns and then its rows.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

temp = rand(size(A,1),size(A,2));
temp1 = sort(temp,1);   % sort along columns
temp2 = sort(temp1,2);  % sort along rows

for c = 1:size(A,2)
    for r = 1:size(A,1)
        [y,x] = find(temp==temp2(r,c));
        rA(r,c) = A(y,x);
    end
end 
