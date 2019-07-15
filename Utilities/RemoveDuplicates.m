function[A]=RemoveDuplicates(A)
% This function removes duplicate row entries from an n x a matrix,
% where 'a' can be any natural number. 

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

for n = 1:size(A,1)-1
    for m = size(A,1):-1:n+1
        if isequal(A(n,:),A(m,:))
            A(m,:)=[];
        end
    end
end