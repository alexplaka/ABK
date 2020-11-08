function rA = RandArray(A)
% This function randomizes the order of an array's contents.
% The array can be either a row or column vector.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

if isempty(find(size(A)==1)) == 0           % if A is a row or column vector
    
    rA = A;             % Initialize and preallocate memory for rA

    if size(A,1) == 1            % Horizontal array
        temp1 = rand(1,size(A,2));
        temp2 = sort(temp1);
        for x = 1:size(A,2)
            rA(x) = A(temp1==temp2(x));
        end 

    elseif size(A,2) == 1       % Vertical array
        temp1 = rand(size(A,1),1);
        temp2 = sort(temp1);
        for x = 1:size(A,1)
            rA(x) = A(temp1==temp2(x));
        end 
    end
    
else
    error('This function randomizes only row or column vectors.');
end
